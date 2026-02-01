module top
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
    input                 i_valid,
    input   [NB_IN-1:0]   i_sample,
    output  [NB_OUT-1:0]  o_sample
);

    //! ShiftRegister model
    reg  signed [NB_IN - 1 : 0]             register    [FIR_LEN-1 : 0]; //! Matrix for registers

    integer ptr1;
    integer ptr2;
    always @(posedge i_clock or negedge i_reset) begin:shiftRegister
        if (!i_reset) begin
        for(ptr1=0;ptr1<FIR_LEN;ptr1=ptr1+1) begin:init
            register[ptr1] <= {(NB_IN){1'b0}};
        end
        end else begin
        if (i_en == 1'b1) begin
            for(ptr2=0;ptr2<FIR_LEN;ptr2=ptr2+1) begin:srmove
                if(ptr2==FIR_LEN-1)
                    if (i_valid)
                        register[ptr2] <= i_sample;
                    else
                        register[ptr2] <= {(NB_IN){1'b0}};
                else
                    register[ptr2] <= register[ptr2+1];
                end   
            end
        end
    end

    wire signed [NB_COEFF - 1 : 0]          coeff       [FIR_LEN-1 : 0];
    // FFE coefficients (S(8,7), floor quantization)
    assign coeff[20] = 8'b11111111; //  -1
    assign coeff[19] = 8'b11111111; //  -1
    assign coeff[18] = 8'b11111111; //  -1
    assign coeff[17] = 8'b00000000; //   0
    assign coeff[16] = 8'b11111111; //  -1
    assign coeff[15] = 8'b00000000; //   0
    assign coeff[14] = 8'b11111111; //  -1
    assign coeff[13] = 8'b00000000; //   0
    assign coeff[12] = 8'b11111111; //  -1
    assign coeff[11] = 8'b00000000; //   0
    assign coeff[10] = 8'b01111111; //  127 (sat)
    assign coeff[9]  = 8'b10000000; // -128 (sat)
    assign coeff[8]  = 8'b11011101; //  -35
    assign coeff[7]  = 8'b11110001; //  -15
    assign coeff[6]  = 8'b11110111; //   -9
    assign coeff[5]  = 8'b11111010; //   -6
    assign coeff[4]  = 8'b00101100; //   44
    assign coeff[3]  = 8'b11110111; //   -9
    assign coeff[2]  = 8'b11111010; //   -6
    assign coeff[1]  = 8'b11111100; //   -4
    assign coeff[0]  = 8'b11111101; //   -3
    
    reg signed [FIR_LEN*NB_IN-1:0]      xk_flat;
    reg signed [FIR_LEN*NB_COEFF-1:0]   coeff_flat;

    integer k;
    always @(*) begin
        for (k = 0; k < FIR_LEN; k = k + 1) begin
            xk_flat[k*NB_IN +: NB_IN] = register[k];
            coeff_flat[k*NB_COEFF +: NB_COEFF] = coeff[k];
        end
    end


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

endmodule;