`timescale 1ns/1ps

module fir_mem_tb;

    // PARAMETERS
    parameter FIR_LEN = 21;
    parameter NB_COEFF = 8;
    parameter NBF_COEFF = 7;
    parameter NB_IN = 18;
    parameter NBF_IN = 15;
    parameter NB_OUT = 18;
    parameter NBF_OUT = 15;


    // CLOCK / CONTROL
    reg clk;  
    reg rst;   
    reg en;
    reg valid;

    // FIR INPUT / OUTPUT
    reg signed  [NB_IN-1:0]     sample_in;
    wire signed [NB_OUT-1:0]    fir_out;

    // FILE HANDLING
    integer fd;
    integer out_fd;
    integer status;
    integer cycle_cnt;

    // CLOCK GENERATION
    always #5 clk = ~clk;   // 100 MHz

    // OPEN INPUT FILE
    initial begin
        fd = $fopen("C:/Users/denis/Documents/beca/porcom-tp-final/python/output_mem/channel_symbols.mem", "r");
        if (fd == 0) begin
            $display("ERROR: cannot open channel_symbols.mem");
            $finish;
        end
    end

    // OPEN OUTPUT FILE
    initial begin
        out_fd = $fopen(
            "C:/Users/denis/Documents/beca/porcom-tp-final/verilog/testbench/fir_out_rtl.txt",
            "w"
        );
        if (out_fd == 0) begin
            $fatal(1, "ERROR: cannot open fir_out_rtl.txt");
        end
    end

    // RESET / ENABLE SEQUENCE
    initial begin
        sample_in = {NB_IN{1'b0}};
        cycle_cnt = 0;

        clk = 1'b0;
        rst = 1'b1;
        en = 1'b0;
        valid = 1'b0;

        #20;
        @(posedge clk);
        rst   = 1'b0;

        #20;
        @(posedge clk);
        rst   = 1'b1;

        #20;
        @(posedge clk);
        en    = 1'b1;
        valid = 1'b1;
    end

    // STREAM INPUT SAMPLES
    // WRITE OUTPUT SAMPLES TO FILE
    always @(posedge clk) begin
        if (!rst) begin
            sample_in <= {NB_IN{1'b0}};
            cycle_cnt <= 0;
        end
        else if (en && valid) begin
            status = $fscanf(fd, "%b\n", sample_in);
            $fwrite(out_fd, "%0d\n", fir_out);
            if (status != 1) begin
                $display("End of channel_symbols.mem at cycle %0d", cycle_cnt);
                $fclose(fd);
                $fclose(out_fd);
                $finish;
            end
            cycle_cnt <= cycle_cnt + 1;
        end
    end


    // DUT INSTANTIATION
    fir #(
        .FIR_LEN(FIR_LEN),
        .NB_COEFF(NB_COEFF),
        .NBF_COEFF(NBF_COEFF),
        .NB_IN(NB_IN),
        .NBF_IN(NBF_IN),
        .NB_OUT(NB_OUT),
        .NBF_OUT(NBF_OUT)
    )
    FIR_DUT (
        .clk       (clk),
        .i_reset   (rst),
        .i_is_data (sample_in),
        .i_en      (en),
        .i_valid   (valid),
        .o_sample  (fir_out)
    );

endmodule
