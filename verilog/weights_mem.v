module weights_mem #(
    parameter FIR_LEN   = 21,
    parameter NB_COEFF  = 8,
    parameter NBF_COEFF = 7,              // Fractional bits
    parameter CENTRAL_TAP = FIR_LEN/2
)(
    input                                   i_clock,
    input                                   i_reset,
    input                                   i_update_en,   // from CMA controller
    input  signed [FIR_LEN*NB_COEFF-1:0]    i_w_new_flat, // CMA-computed weights
    output signed [FIR_LEN*NB_COEFF-1:0]    o_w_flat      // to FIR
);

// Reset value: 1.0 in S(NB_COEFF, NBF_COEFF) = (2^NBF_COEFF - 1) = (2^(NB_COEFF-1) - 1)
localparam RESET_VAL = (1 << (NBF_COEFF)) - 1;  // Represents ~1.0 in the given representation

reg signed [NB_COEFF-1:0] coef [FIR_LEN-1 : 0];
reg signed [FIR_LEN*NB_COEFF-1:0]   coeff_flat;

integer k;
always @(posedge i_clock or negedge i_reset) begin
    if (!i_reset) begin
        for (k = 0; k < FIR_LEN; k = k + 1)
            coef[k] <= {(NB_COEFF){1'b0}};
        // central tap = 1.0 in the fixed-point representation
        coef[CENTRAL_TAP] <= RESET_VAL[NB_COEFF-1:0];
    end
    else if (i_update_en) begin
        for (k = 0; k < FIR_LEN; k = k + 1)
            coef[k] <= i_w_new_flat[k*NB_COEFF +: NB_COEFF];
    end
end

integer j;
always @(*) begin
    for (j = 0; j < FIR_LEN; j = j + 1) begin
        coeff_flat[j*NB_COEFF +: NB_COEFF] = coef[j];
    end
end

assign o_w_flat = coeff_flat;

endmodule