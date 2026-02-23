`timescale 1ns/1ps
// Testbench for parallel top module with PARALLELISM=8
//
// Reads samples from channel_symbols.mem in batches of 8
// Processes 8 samples per clock cycle
// Outputs 8 FIR results per cycle to file

module top_parallelism_tb;

    // ── PARAMETERS ──────────────────────────────────────────────────────────
    parameter FIR_LEN      = 21;
    parameter NB_COEFF     = 28;
    parameter NBF_COEFF    = 23;
    parameter NB_IN        = 18;
    parameter NBF_IN       = 15;
    parameter NB_OUT       = 18;
    parameter NBF_OUT      = 15;
    parameter NB_MU        = 16;
    parameter PARALLELISM  = 8;
    parameter signed [NB_MU-1:0] mu_cma = 16'sd32;

    localparam MEM_LEN = FIR_LEN + PARALLELISM - 1;  // 28

    // ── CLOCK / CONTROL ────────────────────────────────────────────────────
    reg                                     clk;
    reg                                     rst;
    reg                                     en;
    reg                                     valid;

    // ── DUT I/O ────────────────────────────────────────────────────────────
    reg  signed [PARALLELISM*NB_IN-1:0]    sample_in;
    wire signed [PARALLELISM*NB_OUT-1:0]   fir_out;

    // ── UNPACKED VIEWS (for debugging) ─────────────────────────────────────
    reg  signed [NB_IN-1:0]  sample_in_arr  [0:PARALLELISM-1];
    wire signed [NB_OUT-1:0] fir_out_arr    [0:PARALLELISM-1];

    genvar gi;
    generate
        for (gi = 0; gi < PARALLELISM; gi = gi + 1) begin : unpack_fir_o
            assign fir_out_arr[gi] = fir_out[(gi+1)*NB_OUT - 1 -: NB_OUT];
        end
    endgenerate

    integer ui;
    always @(*) begin
        for (ui = 0; ui < PARALLELISM; ui = ui + 1)
            sample_in_arr[ui] = sample_in[(ui+1)*NB_IN - 1 -: NB_IN];
    end

    // ── FILE HANDLING ──────────────────────────────────────────────────────
    integer fd;
    integer out_fd;
    integer coeff_fd;
    integer status;
    integer cycle_cnt;
    integer sample_cnt;

    // ── SAMPLE BUFFER (batch 8 samples at a time) ──────────────────────────
    reg signed [NB_IN-1:0] sample_buffer [0:PARALLELISM-1];
    integer buf_valid;
    integer i, j;

    // ── CLOCK GENERATION (100 MHz) ─────────────────────────────────────────
    always #5 clk = ~clk;   // 100 MHz

    // ── OPEN INPUT FILE ───────────────────────────────────────────────────
    initial begin
        fd = $fopen("C:/Users/denis/Documents/beca/porcom-tp-final/python/output_mem/channel_symbols.mem", "r");
        if (fd == 0) begin
            $display("ERROR: cannot open channel_symbols.mem");
            $finish;
        end
    end

    // ── OPEN OUTPUT FILES ──────────────────────────────────────────────────
    initial begin
        out_fd = $fopen(
            "C:/Users/denis/Documents/beca/porcom-tp-final/verilog/testbench/output/out_ffe_parallel.txt","w");
        coeff_fd = $fopen(
            "C:/Users/denis/Documents/beca/porcom-tp-final/verilog/testbench/output/coeff_ffe_parallel.txt","w");
        if (out_fd == 0) begin
            $fatal(1, "ERROR: cannot open out_ffe_parallel.txt");
        end else if (coeff_fd == 0) begin
            $fatal(1, "ERROR: cannot open coeff_ffe_parallel.txt");
        end
        $fwrite(out_fd, "# Parallel FIR outputs: out[0], out[1], ..., out[7]\n");
        $fwrite(coeff_fd, "# Parallel coefficients: coeff[0], coeff[1], ..., coeff[20]\n");
    end

    // ── RESET / ENABLE SEQUENCE ───────────────────────────────────────────
    initial begin
        sample_in = {PARALLELISM*NB_IN{1'b0}};
        cycle_cnt = 0;
        sample_cnt = 0;

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

    // ── STREAM INPUT SAMPLES (8 per cycle) ─────────────────────────────────
    // WRITE OUTPUT SAMPLES TO FILE
    integer i_coef;
    always @(posedge clk) begin
        if (!rst) begin
            sample_in <= {PARALLELISM*NB_IN{1'b0}};
            cycle_cnt <= 0;
            sample_cnt <= 0;
        end
        else if (en && valid) begin
            // Read up to 8 samples from file
            buf_valid = 0;
            for (i = 0; i < PARALLELISM; i = i + 1) begin
                status = $fscanf(fd, "%b\n", sample_buffer[i]);
                if (status == 1) begin
                    sample_cnt = sample_cnt + 1;
                    buf_valid = i + 1;
                end else begin
                    $display("EOF");
                    $fclose(fd);
                    $fclose(out_fd);
                    $fclose(coeff_fd);
                    $finish;
                end
            end

            // Pack samples into input bus (oldest at LSB, newest at MSB)
            for (j = 0; j < PARALLELISM; j = j + 1) begin
                sample_in[(j+1)*NB_IN - 1 -: NB_IN] <= sample_buffer[j];
            end

            // Write 8 FIR outputs to file
            for (i = 0; i < PARALLELISM; i = i + 1) begin
                $fwrite(out_fd, " %0d", fir_out_arr[i]);
            end
            $fwrite(out_fd, "\n");

            // Write 21 coefficients to file (read from DUT internal memory)
            for (i_coef = 0; i_coef < FIR_LEN; i_coef = i_coef + 1) begin
                $fwrite(coeff_fd, " %0d", top_dut.FFE_MEM.coef[i_coef]);
            end
            $fwrite(coeff_fd, "\n");

            cycle_cnt <= cycle_cnt + 1;

            // Stop when all samples exhausted
            if (buf_valid < PARALLELISM) begin
                $display("End of samples at cycle %0d, total samples read: %0d", cycle_cnt, sample_cnt);
                $fclose(fd);
                $fclose(out_fd);
                $fclose(coeff_fd);
                $finish;
            end
        end
    end

    // ── DUT INSTANTIATION ──────────────────────────────────────────────────
    top #(
        .FIR_LEN       (FIR_LEN),
        .NB_COEFF      (NB_COEFF),
        .NBF_COEFF     (NBF_COEFF),
        .NB_IN         (NB_IN),
        .NBF_IN        (NBF_IN),
        .NB_OUT        (NB_OUT),
        .NBF_OUT       (NBF_OUT),
        .NB_MU         (NB_MU),
        .PARALLELISM   (PARALLELISM)
    )
    top_dut (
        .i_clock  (clk),
        .i_reset  (rst),
        .i_en     (en),
        .i_valid  (valid),
        .i_sample (sample_in),
        .i_mu     (mu_cma),
        .o_sample (fir_out)
    );

endmodule
