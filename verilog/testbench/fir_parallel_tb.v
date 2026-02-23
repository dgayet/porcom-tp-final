`timescale 1ns/1ps
// Testbench for parallel FIR with shift register and weights memory
//
// Reads first 100 samples from channel_symbols_int.mem
// Processes 8 samples per cycle
// Outputs 8 FIR results per cycle to file

module fir_parallel_tb;

    // ── PARAMETERS ──────────────────────────────────────────────────────────
    parameter FIR_LEN      = 21;
    parameter NB_COEFF     = 28;
    parameter NBF_COEFF    = 23;
    parameter NB_IN        = 18;
    parameter NBF_IN       = 15;
    parameter NB_OUT       = 18;
    parameter NBF_OUT      = 15;
    parameter PARALLELISM  = 8;
    parameter CENTRAL_TAP  = 10;

    localparam MEM_LEN = FIR_LEN + PARALLELISM - 1;  // 28

    // ── CLOCK / CONTROL ────────────────────────────────────────────────────
    reg                             clk;
    reg                             rst;
    reg                             en;
    reg                             vld;

    // ── SHIFT REGISTER INPUT / OUTPUT ──────────────────────────────────────
    reg  signed [PARALLELISM*NB_IN-1:0]  shift_i_data;
    wire signed [MEM_LEN*NB_IN-1:0]      shift_o_data;

    // ── WEIGHT MEMORY I/O ──────────────────────────────────────────────────
    wire signed [FIR_LEN*NB_COEFF-1:0]   coeff_flat;
    reg  signed [FIR_LEN*NB_COEFF-1:0]   coeff_in;
    reg                                  coeff_update;

    // ── FIR OUTPUTS ────────────────────────────────────────────────────────
    wire signed [PARALLELISM*NB_OUT-1:0] fir_out_flat;

    // ── UNPACKED VIEWS (for waveform visibility) ───────────────────────────
    reg  signed [NB_IN-1:0]  shift_i_arr  [0:PARALLELISM-1];
    wire signed [NB_IN-1:0]  shift_o_arr  [0:MEM_LEN-1];
    wire signed [NB_OUT-1:0] fir_out_arr  [0:PARALLELISM-1];

    // ── UNPACK SHIFT REGISTER OUTPUT ───────────────────────────────────────
    genvar g0;
    generate
        for (g0 = 0; g0 < MEM_LEN; g0 = g0 + 1) begin : unpack_shift_o
            assign shift_o_arr[g0] = shift_o_data[(g0+1)*NB_IN - 1 -: NB_IN];
        end
    endgenerate

    // ── UNPACK FIR OUTPUT ──────────────────────────────────────────────────
    genvar g1;
    generate
        for (g1 = 0; g1 < PARALLELISM; g1 = g1 + 1) begin : unpack_fir_o
            assign fir_out_arr[g1] = fir_out_flat[(g1+1)*NB_OUT - 1 -: NB_OUT];
        end
    endgenerate

    // ── COMBINATIONAL UNPACK OF SHIFT INPUT ────────────────────────────────
    integer u0;
    always @(*) begin
        for (u0 = 0; u0 < PARALLELISM; u0 = u0 + 1)
            shift_i_arr[u0] = shift_i_data[(u0+1)*NB_IN - 1 -: NB_IN];
    end

    // ── SHIFT REGISTER INSTANCE ────────────────────────────────────────────
    shift_register #(
        .N          (FIR_LEN),
        .NB         (NB_IN),
        .PARALLELISM(PARALLELISM)
    ) shift_reg_inst (
        .i_clock (clk),
        .i_reset (rst),
        .i_enable(en),
        .i_valid (vld),
        .i_data  (shift_i_data),
        .o_data  (shift_o_data)
    );

    // ── WEIGHTS MEMORY INSTANCE ────────────────────────────────────────────
    weights_mem #(
        .FIR_LEN    (FIR_LEN),
        .NB_COEFF   (NB_COEFF),
        .NBF_COEFF  (NBF_COEFF),
        .CENTRAL_TAP(CENTRAL_TAP)
    ) weights_mem_inst (
        .i_clock    (clk),
        .i_reset    (rst),
        .i_update_en(coeff_update),
        .i_w_new_flat(coeff_in),
        .o_w_flat   (coeff_flat)
    );

    // ── 8 PARALLEL FIR INSTANCES ───────────────────────────────────────────
    genvar fir_i;
    generate
        for (fir_i = 0; fir_i < PARALLELISM; fir_i = fir_i + 1) begin : fir_gen
            fir #(
                .FIR_LEN   (FIR_LEN),
                .NB_COEFF  (NB_COEFF),
                .NBF_COEFF (NBF_COEFF),
                .NB_IN     (NB_IN),
                .NBF_IN    (NBF_IN),
                .NB_OUT    (NB_OUT),
                .NBF_OUT   (NBF_OUT)
            ) fir_inst (
                .o_sample  (fir_out_flat[(fir_i+1)*NB_OUT - 1 -: NB_OUT]),
                .clk       (clk),
                .i_reset   (rst),
                .i_data_reg(shift_o_data[(fir_i+FIR_LEN)*NB_IN - 1 -: FIR_LEN*NB_IN]),
                .i_coeff   (coeff_flat),
                .i_en      (en),
                .i_valid   (vld)
            );
        end
    endgenerate

    // ── CLOCK GENERATION (10 ns period = 100 MHz) ──────────────────────────
    always #5 clk = ~clk;

    // ── FILE HANDLES ───────────────────────────────────────────────────────
    integer fd_in;
    integer status;
    integer sample_cnt;
    integer cycle_cnt;

    // ── SAMPLE BUFFER FOR BATCH-MODE ──────────────────────────────────────
    reg signed [NB_IN-1:0] sample_buffer [0:7];  // up to 8 samples per cycle

    // ── OPEN INPUT FILE ────────────────────────────────────────────────────
    initial begin
        fd_in = $fopen("C:/Users/denis/Documents/beca/porcom-tp-final/python/output_mem/channel_symbols_int.mem", "r");
        if (fd_in == 0) begin
            $display("ERROR: cannot open channel_symbols_int.mem");
            $finish;
        end
    end

    // ── RESET / ENABLE SEQUENCE ───────────────────────────────────────────
    initial begin
        clk = 1'b0;
        rst = 1'b1;
        en  = 1'b0;
        vld = 1'b0;
        coeff_update = 1'b0;
        coeff_in = {FIR_LEN*NB_COEFF{1'b0}};
        shift_i_data = {PARALLELISM*NB_IN{1'b0}};

        #20;
        @(posedge clk);
        rst = 1'b0;

        #20;
        @(posedge clk);
        rst = 1'b1;

        #20;
        @(posedge clk);
        en  = 1'b1;
        vld = 1'b1;
    end

    // ── MAIN TEST LOOP ────────────────────────────────────────────────────
    integer i, j;
    always @(posedge clk) begin
        if (!rst) begin
            shift_i_data <= {PARALLELISM*NB_IN{1'b0}};
            sample_cnt <= 0;
            cycle_cnt <= 0;
        end
        else if (en && vld) begin
            // Try to fill buffer up to 8 samples (or stop at 100)
            for (i = 0; i < PARALLELISM; i = i + 1) begin
                if (sample_cnt < 100) begin
                    status = $fscanf(fd_in, "%d\n", sample_buffer[i]);
                    if (status == 1) begin
                        sample_cnt = sample_cnt + 1;
                    end else begin
                        $display("Reached EOF");
                        $fclose(fd_in);
                        $finish;
                    end
                end
            end

            // Pack the batch (in the order read: oldest=first)
            for (j = 0; j < PARALLELISM; j = j + 1) begin
                shift_i_data[(j+1)*NB_IN - 1 -: NB_IN] <= sample_buffer[j];
            end

            cycle_cnt = cycle_cnt + 1;

            // Stop after 100 samples have been read
            if (sample_cnt >= 100) begin
                $display("Processed 100 samples in %0d cycles", cycle_cnt);
                #20;
                $fclose(fd_in);
                $finish;
            end
        end
    end

endmodule
