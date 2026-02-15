// CMA Full Module - Updates all FFE taps
// Instantiates cma_base for each tap and manages error computation + startup delay
// Convention:
// i_coeff[k] and i_xk[k] are indexed where k=FFE_LEN-1 is newest (x[N-1], h[N-1])
// i_coeff[0] and i_xk[0] are oldest samples
// Matches Python CMA implementation

module cma_full
#(
    parameter NB_I      = 18,
    parameter NBF_I     = 15,
    parameter FFE_LEN   = 21,
    parameter NB        = 8,
    parameter NBF       = 7,
    parameter NBF_COEFF = 23,    // Fractional bits of cma_r input
    parameter NB_MU     = 16
)
(
    input                               i_clock,
    input                               i_reset,
    input                               i_valid,
    input signed [NB_I-1:0]             i_fir_out,
    input signed [NB-1:0]               cma_r,        // CMA radius (e.g., 0.5125)
    input signed [NB_I*FFE_LEN-1:0]     i_xk_flat,
    input signed [FFE_LEN*NB-1:0]       i_coeff_flat,
    input signed [NB_MU-1:0]            i_mu,
    output signed [FFE_LEN*NB-1:0]      o_new_coeff
);

    localparam STARTUP_DELAY = 3*FFE_LEN;

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
    // CMA Error Computation
    // ============================================================================
    
    wire signed [2*NB_I-1:0] y_sq_full;
    assign y_sq_full = i_fir_out * i_fir_out;

    // Scale CMA_R by the difference in fractional bits
    // CRITICAL: Sign-extend cma_r to full width BEFORE shifting to avoid overflow
    // y_sq_full has 2*NBF_I fractional bits, cma_r has NBF_COEFF fractional bits
    // Need to shift left by (2*NBF_I - NBF_COEFF) to align them
    wire signed [2*NB_I-1:0] R_full;
    assign R_full = {{(2*NB_I - NB){cma_r[NB-1]}}, cma_r} <<< (2*NBF_I - NBF);

    // Error = y^2 - R, scaled back to coefficient fractional format
    wire signed [2*NB_I:0] err_full;
    wire signed [NB-1:0] error;
    assign err_full = y_sq_full - R_full;
    assign error = err_full >>> (2*NBF_I - NBF);

    // ============================================================================
    // Instantiate CMA Base Updaters for each tap
    // ============================================================================
    
    genvar tap_idx;
    generate
        for (tap_idx = 0; tap_idx < FFE_LEN; tap_idx = tap_idx + 1) begin : cma_tap_gen
            cma_base #(
                .NB_I(NB_I),
                .NBF_I(NBF_I),
                .NB(NB),
                .NBF(NBF),
                .NB_MU(NB_MU)
            ) cma_tap (
                .i_xk(xk[tap_idx]),
                .i_fir_out(i_fir_out),
                .i_error(error),
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

    assign o_new_coeff = w_new_flat; // Hold old coeffs until update enabled

endmodule
