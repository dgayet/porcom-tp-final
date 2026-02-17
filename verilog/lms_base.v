// LMS Base Module - Updates a single tap
// Computes: w_new[k] = w[k] - mu * error * x[k]
// Error is provided at full precision (NB_I+1 bits with NBF_I fractional bits)

module lms_base
#(
    parameter NB_I      = 18,      // x[k], error width (input width)
    parameter NBF_I     = 15,      // x[k], error fractional bits
    parameter NB_ERROR  = 19,      // error width (NB_I+1 for subtraction result)
    parameter NBF_ERROR = 15,      // error fractional bits (same as NBF_I)
    parameter NB        = 8,       // coefficient width
    parameter NBF       = 7,       // coefficient fractional bits
    parameter NB_MU     = 16       // mu width (assumed signed, frac = NB_MU-1)
)
(
    input  signed [NB_I-1:0]         i_xk,
    input  signed [NB_ERROR-1:0]     i_error,
    input  signed [NB-1:0]           i_w,
    input  signed [NB_MU-1:0]        i_mu,
    output signed [NB-1:0]           o_w_new
);

    // -----------------------------
    // Stage 1: error * mu
    // error: NB_ERROR bits, NBF_ERROR fractional bits
    // mu: NB_MU bits, NB_MU-1 fractional bits
    // Result: NB_ERROR + NB_MU bits, NBF_ERROR + NB_MU-1 fractional bits
    // -----------------------------
    localparam NB_M1  = NB_ERROR + NB_MU;
    localparam NBF_M1 = NBF_ERROR + (NB_MU-1);

    wire signed [NB_M1-1:0] m1;
    assign m1 = i_error * i_mu;

    // -----------------------------
    // Stage 2: (error * mu) * xk
    // m1: NB_M1 bits, NBF_M1 fractional bits
    // xk: NB_I bits, NBF_I fractional bits
    // Result: NB_M1 + NB_I bits, NBF_M1 + NBF_I fractional bits
    // -----------------------------
    localparam NB_UPD_FULL  = NB_M1 + NB_I;
    localparam NBF_UPD_FULL = NBF_M1 + NBF_I;

    wire signed [NB_UPD_FULL-1:0] upd_full;
    assign upd_full = m1 * i_xk;

    // -----------------------------
    // Extend w to full precision
    // w: NB bits, NBF fractional bits
    // Need to shift to NBF_UPD_FULL and extend to NB_UPD_FULL
    // Left shift by (NBF_UPD_FULL - NBF) positions
    // -----------------------------
    wire signed [NB_UPD_FULL-1:0] w_ext;
    assign w_ext =
        {{(NB_UPD_FULL-NB){i_w[NB-1]}}, i_w}
        <<< (NBF_UPD_FULL - NBF);

    // -----------------------------
    // Full precision subtraction
    // w_ext and upd_full both have NBF_UPD_FULL fractional bits
    // -----------------------------
    wire signed [NB_UPD_FULL-1:0] w_new_full;
    assign w_new_full = w_ext - upd_full;

    // -----------------------------
    // Saturation + truncation
    // Extract back to NB-bit coefficient with NBF fractional bits
    // Need to shift right by (NBF_UPD_FULL - NBF) positions
    // Check for overflow in integer bits
    // -----------------------------
    localparam NBI_FULL   = NB_UPD_FULL - NBF_UPD_FULL;
    localparam NBI_OUT    = NB - NBF;
    localparam NB_SAT     = NBI_FULL - NBI_OUT;

    wire of_pos = ~|w_new_full[NB_UPD_FULL-1 -: NB_SAT+1];
    wire of_neg =  &w_new_full[NB_UPD_FULL-1 -: NB_SAT+1];

    wire signed [NB-1:0] w_new_sat =
        (of_pos || of_neg) ?
            w_new_full[NB_UPD_FULL-(NBI_FULL-NBI_OUT)-1 -: NB] :
        (w_new_full[NB_UPD_FULL-1]) ?
            {1'b1,{NB-1{1'b0}}} :
            {1'b0,{NB-1{1'b1}}};

    assign o_w_new = w_new_sat;

endmodule
