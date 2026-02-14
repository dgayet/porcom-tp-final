// Convention:
// coeff[0] = newest tap (x[n])
// coeff[k] multiplies x[n-k]
// Matches Python CMA implementation

module top_old
#(
    parameter FIR_LEN   = 21,
    parameter NB_COEFF  = 8,
    parameter NBF_COEFF = 7,
    parameter NB_IN     = 8,
    parameter NBF_IN    = 7,
    parameter NB_OUT    = 18,
    parameter NBF_OUT   = 17
)
(
    input                 i_clock,
    input                 i_reset,
    input                 i_en,
    input                 i_update_en,
    input                 i_valid,
    input   [NB_IN-1:0]   i_sample,
    output  [NB_OUT-1:0]  o_sample
);

    // THIS WILL BE REPLACED WITH CMA
    wire signed [NB_COEFF - 1 : 0]          coeff       [FIR_LEN-1 : 0];
    // FFE coefficients (S(8,7), floor quantization)
    assign coeff[0] = 8'b11111111; //  -1
    assign coeff[1] = 8'b11111111; //  -1
    assign coeff[2] = 8'b11111111; //  -1
    assign coeff[3] = 8'b00000000; //   0
    assign coeff[4] = 8'b11111111; //  -1
    assign coeff[5] = 8'b00000000; //   0
    assign coeff[6] = 8'b11111111; //  -1
    assign coeff[7] = 8'b00000000; //   0
    assign coeff[8] = 8'b11111111; //  -1
    assign coeff[9] = 8'b00000000; //   0
    assign coeff[10] = 8'b01111111; //  127 (sat)
    assign coeff[11]  = 8'b10000000; // -128 (sat)
    assign coeff[12]  = 8'b11011101; //  -35
    assign coeff[13]  = 8'b11110001; //  -15
    assign coeff[14]  = 8'b11110111; //   -9
    assign coeff[15]  = 8'b11111010; //   -6
    assign coeff[16]  = 8'b00101100; //   44
    assign coeff[17]  = 8'b11110111; //   -9
    assign coeff[18]  = 8'b11111010; //   -6
    assign coeff[19]  = 8'b11111100; //   -4
    assign coeff[20]  = 8'b11111101; //   -3
    
    reg signed [FIR_LEN*NB_COEFF-1:0]   new_coeff_flat;
    integer j;
    always @(*) begin
        for (j = 0; j < FIR_LEN; j = j + 1) begin
            new_coeff_flat[j*NB_COEFF +: NB_COEFF] = coeff[j];
        end
    end
    
    //! Weights registering
    wire signed [FIR_LEN*NB_COEFF-1:0]   coeff_flat;
    
    weights_mem #(
        .FIR_LEN(FIR_LEN),
        .NB_COEFF(NB_COEFF),
        .CENTRAL_TAP(10)
    ) FFE_MEM (
        .i_clock(i_clock),
        .i_reset(i_reset),
        .i_update_en(i_update_en),
        .i_w_new_flat(new_coeff_flat),
        .o_w_flat(coeff_flat)
    );

    //! ShiftRegister model
    wire signed [FIR_LEN*NB_IN-1:0]  xk_flat;
    shift_register #(
        .N(FIR_LEN),
        .NB(NB_IN)
    ) DELAY_LINE (
        .i_clock(i_clock),
        .i_reset(i_reset),
        .i_enable(i_en),
        .i_valid(i_valid),
        .i_data(i_sample),
        .o_data(xk_flat)
    );

    fir #(
        .FIR_LEN(FIR_LEN),
        .NB_COEFF(NB_COEFF),
        .NBF_COEFF(NBF_COEFF),
        .NB_IN(NB_IN),
        .NBF_IN(NBF_IN),
        .NB_OUT(NB_OUT),
        .NBF_OUT(NBF_OUT)
    ) fir_inst (
        .o_sample(o_sample),
        .clk(i_clock),
        .i_reset(i_reset),
        .i_data_reg(xk_flat),
        .i_coeff(coeff_flat),
        .i_en(i_en),
        .i_valid(i_valid)
    );

endmodule