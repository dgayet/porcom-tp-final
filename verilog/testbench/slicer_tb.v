`timescale 1ns/100ps

module slicer_tb();
    parameter NB = 8;
    parameter NF = 7;
    
    reg                     clk;
    reg                     reset;
    reg                     en;
    reg                     valid;
    reg  signed [NB-1:0]    sample_in;
    wire signed [NB-1:0]    slicer_out;
    wire        [1:0]       gray_level;


    initial begin
        clk = 1'b0;
        reset = 1'b1;
        en = 1'b0;
        valid = 1'b0;
        sample_in = 8'b01111111;
        #100;
        @(posedge clk);
        reset = 1'b0;
        #100;
        reset = 1'b1;
        #100;
        en = 1'b1;
        valid = 1'b1;
        #100;
        @(posedge clk);
        sample_in = 8'b11111111;
        #100;
        @(posedge clk);
        sample_in = 8'b10000000;
        #100;
        @(posedge clk);
        sample_in = 8'b00011100;
        #100;    
        @(posedge clk);
        sample_in = 8'b11110000;
        #100;     
        reset = 1'b0;
        #100;   
    

        $finish;
    end

    // Clock generation
    always #5 clk = ~clk;

    // Instantiate the slicer_pam4 module
    slicer_pam4 #(
        .NB(NB),
        .NF(NF)
    ) DUT (
        .i_clock(clk),
        .i_reset(reset),
        .i_enable(en),
        .i_valid(valid),
        .i_sample(sample_in),
        .o_slicer(slicer_out),
        .o_gray_level(gray_level)
    );


endmodule