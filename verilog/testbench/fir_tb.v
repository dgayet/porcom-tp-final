`timescale 1ns/100ps

module fir_tb();
    parameter FIR_LEN = 24;
    parameter NB_COEFF = 8;
    parameter NBF_COEFF = 7;
    parameter NB_IN = 18;
    parameter NBF_IN = 15;
    parameter NB_OUT = 18;
    parameter NBF_OUT = 15;

    reg clk;
    reg reset;
    reg en;
    reg valid;

    reg signed  [NB_IN-1:0]     input_sample;
    wire signed [NB_OUT-1:0]    output_sample;

    initial begin
        clk = 1'b0;
        reset = 1'b1;
        en = 1'b0;
        valid = 1'b0;
        input_sample = 18'b001000000000000000;
        #100;
        @(posedge clk);
        reset = 1'b0;
        #100;
        reset = 1'b1;
        #100;
        en = 1'b1;
        valid = 1'b1;
        #5;
        @(posedge clk);
        input_sample = 18'b000000000000000000;
        #120;
        @(posedge clk);
        input_sample = 18'b001000000000000000;
        #200;
        @(posedge clk);
        reset = 1'b0;
        #100;
        
        $finish;
    end

    // clk generation
    always #5 clk = ~clk;

    // module instantiation
    fir 
    #(
        .FIR_LEN(FIR_LEN),
        .NB_COEFF(NB_COEFF),
        .NBF_COEFF(NBF_COEFF),
        .NB_IN(NB_IN),
        .NBF_IN(NBF_IN),
        .NB_OUT(NB_OUT),
        .NBF_OUT(NBF_OUT)
    )
    FIR_DUT (
        .o_sample(output_sample),
        .clk(clk),
        .i_reset(reset),
        .i_is_data(input_sample),
        .i_en(en),
        .i_valid(valid)
    );

endmodule
