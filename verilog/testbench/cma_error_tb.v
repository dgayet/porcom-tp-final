`timescale 1ns/1ps

module cma_error_tb;
    // Latency summary:
    // input sample x[n]        @ cycle n
    // FIR output y[n]          @ cycle n+1
    // slicer output            @ cycle n+2

    // PARAMETERS
    parameter FIR_LEN = 21;
    parameter NB_COEFF = 8;
    parameter NBF_COEFF = 7;
    parameter NB_IN = 18;
    parameter NBF_IN = 15;
    parameter NB_OUT = 18;
    parameter NBF_OUT = 15;
    parameter signed [NB_COEFF-1:0]   CMA_R = 8'sd65;  // 0.5125 
    parameter NB_MU = 16;


    // CLOCK / CONTROL
    reg clk;  
    reg rst;   
    reg en;
    reg update_en;
    reg valid;

    // FIR INPUT / OUTPUT
    reg signed  [NB_IN-1:0]     sample_in;
    wire signed [NB_OUT-1:0]    fir_out;

    // CMA OUTPUT
    wire signed [NB_COEFF-1:0]    error_out;

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
            "C:/Users/denis/Documents/beca/porcom-tp-final/verilog/testbench/cma_error_rtl.txt",
            "w"
        );
        if (out_fd == 0) begin
            $fatal(1, "ERROR: cannot open cma_error_rtl.txt");
        end
    end

    // RESET / ENABLE SEQUENCE
    initial begin
        sample_in = {NB_IN{1'b0}};
        cycle_cnt = 0;

        clk = 1'b0;
        rst = 1'b1;
        en = 1'b0;
        update_en = 1'b0;
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
        update_en = 1'b1;
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
            $fwrite(out_fd, "%0d\n", error_out);
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
    top_old #(
        .FIR_LEN(FIR_LEN),
        .NB_COEFF(NB_COEFF),
        .NBF_COEFF(NBF_COEFF),
        .NB_IN(NB_IN),
        .NBF_IN(NBF_IN),
        .NB_OUT(NB_OUT),
        .NBF_OUT(NBF_OUT)
    )
    FIR_DUT (
        .i_clock   (clk),
        .i_reset   (rst),
        .i_en      (en),
        .i_update_en(update_en),
        .i_valid   (valid),
        .i_sample (sample_in),
        .o_sample  (fir_out)
    );

    // Instantiate the cma_update_module module
    cma_update #(
        .NB_I(NB_OUT),
        .NBF_I(NBF_OUT),
        .FFE_LEN(FIR_LEN),
        .NB(NB_COEFF),
        .NBF(NBF_COEFF),
        .NB_MU(NB_MU)
    ) DUT (
        .i_clock(clk),
        .i_reset(rst),
        .i_fir_out(fir_out),
        .cma_r(CMA_R),
        .cma_error(error_out)
    );

endmodule
