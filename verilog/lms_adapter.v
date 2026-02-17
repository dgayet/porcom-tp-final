// LMS Adapter Module - Updates all FFE taps
// Instantiates lms_base for each tap and manages coefficient updates
// Convention:
// i_coeff[k] and i_xk[k] are indexed where k=FFE_LEN-1 is newest (x[N-1], h[N-1])
// i_coeff[0] and i_xk[0] are oldest samples
// Matches Python LMS implementation

module lms_adapter
#(
    parameter NB_I      = 18,   // x[k] width
    parameter NBF_I     = 15,   // x[k] fractional bits
    parameter FFE_LEN   = 21,   // Number of FFE taps
    parameter NB        = 8,    // coefficient width
    parameter NBF       = 7,    // coefficient fractional bits
    parameter NB_MU     = 16    // mu width
)
(
    input                               i_clock,
    input                               i_reset,
    input                               i_valid,
    input signed [NB_I-1:0]             i_fir_out,
    input signed [NB_I-1:0]             i_slicer_out,
    input signed [NB_I*FFE_LEN-1:0]     i_xk_flat,
    input signed [FFE_LEN*NB-1:0]       i_coeff_flat,
    input signed [NB_MU-1:0]            i_mu,
    output signed [FFE_LEN*NB-1:0]      o_new_coeff
);

    // Unpacking arrays
    reg signed [NB_I-1:0] xk [0:FFE_LEN-1];
    reg signed [NB-1:0]  w  [0:FFE_LEN-1];
    wire signed [NB-1:0] w_new [0:FFE_LEN-1];
    reg signed [FFE_LEN*NB-1:0] w_new_flat;

    integer k;

    // Unpack flat signals
    always @(*) begin
        for (k = 0; k < FFE_LEN; k = k + 1) begin
            xk[k] = i_xk_flat[k*NB_I +: NB_I];
            w[k]  = i_coeff_flat[k*NB +: NB];
        end
    end

    // ============================================================================
    // LMS Error Computation
    // ============================================================================
    // Error has full precision: NB_I+1 bits (for subtraction), NBF_I fractional bits
    
    wire signed [NB_I:0] error_full;
    assign error_full = i_fir_out - i_slicer_out;

    // ============================================================================
    // Instantiate LMS Base Updaters for each tap
    // ============================================================================
    
    genvar tap_idx;
    generate
        for (tap_idx = 0; tap_idx < FFE_LEN; tap_idx = tap_idx + 1) begin : lms_tap_gen
            lms_base #(
                .NB_I(NB_I),
                .NBF_I(NBF_I),
                .NB_ERROR(NB_I+1),
                .NBF_ERROR(NBF_I),
                .NB(NB),
                .NBF(NBF),
                .NB_MU(NB_MU)
            ) lms_tap (
                .i_xk(xk[tap_idx]),
                .i_error(error_full),
                .i_w(w[tap_idx]),
                .i_mu(i_mu),
                .o_w_new(w_new[tap_idx])
            );
        end
    endgenerate

    // ============================================================================
    // Flatten updated coefficients
    // ============================================================================
    
    always @(*) begin
        for (k = 0; k < FFE_LEN; k = k + 1)
            w_new_flat[k*NB +: NB] = w_new[k];
    end

    assign o_new_coeff = w_new_flat;

endmodule
