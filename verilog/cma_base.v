// CMA Base Module - Updates a single tap
// Computes: w_new[k] = w[k] - mu * error * y * x[k]

module cma_base
#(
    parameter NB_I   = 18,      // x[k], y width
    parameter NBF_I  = 15,      // x[k], y fractional bits
    parameter NB     = 8,       // coefficient width
    parameter NBF    = 7,       // coefficient fractional bits
    parameter NB_MU  = 16       // mu width (assumed signed, frac = NB_MU-1)
)
(
    input  signed [NB_I-1:0]   i_xk,
    input  signed [NB_I-1:0]   i_fir_out,
    input  signed [NB-1:0]     i_error,
    input  signed [NB-1:0]     i_w,
    input  signed [NB_MU-1:0]  i_mu,
    output signed [NB-1:0]     o_w_new
);

    // -----------------------------
    // Stage 1: y * error
    // -----------------------------
    localparam NB_M1  = NB_I + NB;
    localparam NBF_M1 = NBF_I + NBF;

    wire signed [NB_M1-1:0] m1;
    assign m1 = i_fir_out * i_error;

    // -----------------------------
    // Stage 2: (y * error) * mu
    // -----------------------------
    localparam NB_M2  = NB_M1 + NB_MU;
    localparam NBF_M2 = NBF_M1 + (NB_MU-1);

    wire signed [NB_M2-1:0] m2;
    assign m2 = m1 * i_mu;

    // -----------------------------
    // Stage 3: (y * error * mu) * xk
    // -----------------------------
    localparam NB_UPD_FULL  = NB_M2 + NB_I;
    localparam NBF_UPD_FULL = NBF_M2 + NBF_I;

    wire signed [NB_UPD_FULL-1:0] upd_full;
    assign upd_full = m2 * i_xk;

    // -----------------------------
    // Extend w to full precision
    // -----------------------------
    wire signed [NB_UPD_FULL-1:0] w_ext;
    assign w_ext =
        {{(NB_UPD_FULL-NB){i_w[NB-1]}}, i_w}
        <<< (NBF_UPD_FULL - NBF);

    // -----------------------------
    // Full precision subtraction
    // -----------------------------
    wire signed [NB_UPD_FULL-1:0] w_new_full;
    assign w_new_full = w_ext - upd_full;

    // -----------------------------
    // Saturation + truncation
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
