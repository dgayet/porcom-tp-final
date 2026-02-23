`timescale 1ns/1ps
// Testbench for shift_register with PARALLELISM=8

module shift_reg_parallel_tb;

    // ── Parameters ─────────────────────────────────────────────────────────
    localparam N       = 21;            // FIR tap count
    localparam NB      = 18;            // word width
    localparam P       = 8;             // parallelism
    localparam MEM_LEN = N + P - 1;     // 28 slots stored

    // ── DUT Signals ────────────────────────────────────────────────────────
    reg                     clk, rst_n, en, vld;
    reg  [P*NB-1:0]         i_data;
    wire [MEM_LEN*NB-1:0]   o_data;

    // ── Unpacked views (visible as arrays in simulator waveform windows) ────
    // i_data_arr[0] = oldest sample of current batch
    // i_data_arr[P-1] = newest sample of current batch
    reg signed [NB-1:0] i_data_arr [0:P-1];
    // o_data_arr[0] = globally oldest sample in memory
    // o_data_arr[MEM_LEN-1] = globally newest sample in memory
    wire signed [NB-1:0] o_data_arr [0:MEM_LEN-1];

    genvar gi;
    generate
        for (gi = 0; gi < MEM_LEN; gi = gi + 1) begin : unpack_o
            assign o_data_arr[gi] = o_data[gi*NB +: NB];
        end
    endgenerate

    integer ui;
    always @(*) begin
        for (ui = 0; ui < P; ui = ui + 1)
            i_data_arr[ui] = i_data[ui*NB +: NB];
    end

    // ── DUT ────────────────────────────────────────────────────────────────
    shift_register #(
        .N          (N),
        .NB         (NB),
        .PARALLELISM(P)
    ) dut (
        .i_clock(clk),
        .i_reset(rst_n),
        .i_enable(en),
        .i_valid(vld),
        .i_data (i_data),
        .o_data (o_data)
    );

    // ── Clock: 10 ns period ─────────────────────────────────────────────────
    initial clk = 0;
    always #5 clk = ~clk;

    // ── Helper: read one slot as a signed integer ────────────────────────────
    function signed [NB-1:0] slot;
        input integer idx;
        slot = o_data[idx*NB +: NB];
    endfunction

    // ── Task: drive one batch and advance one clock ──────────────────────────
    // start_val: value of the oldest sample in the batch
    // i_data[j*NB +: NB] = start_val + j  (j=0 oldest, j=P-1 newest)
    task apply_batch;
        input integer start_val;
        integer j;
        begin
            for (j = 0; j < P; j = j + 1)
                i_data[j*NB +: NB] = start_val[NB-1:0] + j[NB-1:0];
            vld = 1;
            @(posedge clk); #1;
            vld = 0;
        end
    endtask

    // ── Main test ────────────────────────────────────────────────────────────
    integer k;
    integer errors;

    initial begin
        errors = 0;
        en = 1; vld = 0; i_data = 0;

        // ────────────────────────────────────────────────────────────────────
        // 1. RESET
        // ────────────────────────────────────────────────────────────────────
        rst_n = 0;
        @(posedge clk); #1;
        rst_n = 1;

        for (k = 0; k < MEM_LEN; k = k + 1) begin
            if (slot(k) !== 0) begin
                $error("[RESET] slot[%0d] = %0d, expected 0", k, slot(k));
                errors = errors + 1;
            end
        end
        $display("[RESET] done.");

        // ────────────────────────────────────────────────────────────────────
        // 2. GATE CHECK: memory must not change when valid or enable is low
        // ────────────────────────────────────────────────────────────────────
        i_data = {P*NB{1'b1}};

        // enable=1, valid=0 → no shift
        en = 1; vld = 0;
        @(posedge clk); #1;
        for (k = 0; k < MEM_LEN; k = k + 1) begin
            if (slot(k) !== 0) begin
                $error("[GATE valid=0] slot[%0d] changed", k);
                errors = errors + 1;
            end
        end

        // enable=0, valid=1 → no shift
        en = 0; vld = 1;
        @(posedge clk); #1;
        for (k = 0; k < MEM_LEN; k = k + 1) begin
            if (slot(k) !== 0) begin
                $error("[GATE enable=0] slot[%0d] changed", k);
                errors = errors + 1;
            end
        end

        en = 1; vld =0; i_data = 0;
        $display("[GATE] done.");

        // ────────────────────────────────────────────────────────────────────
        // 3. CYCLE 1 – batch [1..8]
        //    Expected: slots[0..19]=0, slots[20]=1 .. slots[27]=8
        // ────────────────────────────────────────────────────────────────────
        apply_batch(1);

        for (k = 0; k < MEM_LEN - P; k = k + 1) begin
            if (slot(k) !== 0) begin
                $error("[CY1] slot[%0d] = %0d, expected 0", k, slot(k));
                errors = errors + 1;
            end
        end
        for (k = 0; k < P; k = k + 1) begin
            if (slot(MEM_LEN - P + k) !== k + 1) begin
                $error("[CY1] slot[%0d] = %0d, expected %0d",
                       MEM_LEN-P+k, slot(MEM_LEN-P+k), k+1);
                errors = errors + 1;
            end
        end
        $display("[CYCLE 1] oldest-slot(20)=%0d newest-slot(27)=%0d",
                 slot(MEM_LEN-P), slot(MEM_LEN-1));

        // ────────────────────────────────────────────────────────────────────
        // 4. CYCLE 2 – batch [9..16]
        //    Expected: slots[0..11]=0, slots[12..19]=1..8, slots[20..27]=9..16
        // ────────────────────────────────────────────────────────────────────
        apply_batch(9);

        for (k = 0; k < MEM_LEN - 2*P; k = k + 1) begin
            if (slot(k) !== 0) begin
                $error("[CY2] slot[%0d] = %0d, expected 0", k, slot(k));
                errors = errors + 1;
            end
        end
        for (k = 0; k < P; k = k + 1) begin
            if (slot(MEM_LEN - 2*P + k) !== k + 1) begin
                $error("[CY2] old-batch slot[%0d] = %0d, expected %0d",
                       MEM_LEN-2*P+k, slot(MEM_LEN-2*P+k), k+1);
                errors = errors + 1;
            end
        end
        for (k = 0; k < P; k = k + 1) begin
            if (slot(MEM_LEN - P + k) !== k + 9) begin
                $error("[CY2] new-batch slot[%0d] = %0d, expected %0d",
                       MEM_LEN-P+k, slot(MEM_LEN-P+k), k+9);
                errors = errors + 1;
            end
        end
        $display("[CYCLE 2] slot(12)=%0d  slot(19)=%0d  slot(20)=%0d  slot(27)=%0d",
                 slot(MEM_LEN-2*P), slot(MEM_LEN-P-1),
                 slot(MEM_LEN-P),   slot(MEM_LEN-1));

        // ────────────────────────────────────────────────────────────────────
        // 5. CYCLE 3 – batch [17..24]
        // ────────────────────────────────────────────────────────────────────
        apply_batch(17);
        $display("[CYCLE 3] slot(4)=%0d  slot(11)=%0d  slot(20)=%0d  slot(27)=%0d",
                 slot(MEM_LEN-3*P), slot(MEM_LEN-2*P-1),
                 slot(MEM_LEN-P),   slot(MEM_LEN-1));

        // ────────────────────────────────────────────────────────────────────
        // 6. CYCLE 4 – batch [25..32]
        //    Window fully populated. 32 samples pushed, MEM_LEN=28 retained.
        //    slot[k] = k + 5  for k = 0..27
        // ────────────────────────────────────────────────────────────────────
        apply_batch(25);

        for (k = 0; k < MEM_LEN; k = k + 1) begin
            if (slot(k) !== k + 5) begin
                $error("[CY4] slot[%0d] = %0d, expected %0d",
                       k, slot(k), k+5);
                errors = errors + 1;
            end
        end
        $display("[CYCLE 4] slot(0)=%0d (exp 5)  slot(27)=%0d (exp 32)",
                 slot(0), slot(MEM_LEN-1));

        // ────────────────────────────────────────────────────────────────────
        // 7. FIR WINDOW SPOT-CHECK
        //    FIR[i=0] oldest output: uses slots[0..20]  → samples 5..25
        //    FIR[i=7] newest output: uses slots[7..27]  → samples 12..32
        // ────────────────────────────────────────────────────────────────────
        if (slot(7)  !== 12) begin
            $error("[FIR window] FIR[7] low  slot[7]=%0d, expected 12",  slot(7));
            errors = errors + 1;
        end
        if (slot(27) !== 32) begin
            $error("[FIR window] FIR[7] high slot[27]=%0d, expected 32", slot(27));
            errors = errors + 1;
        end
        if (slot(0)  !== 5) begin
            $error("[FIR window] FIR[0] low  slot[0]=%0d, expected 5",   slot(0));
            errors = errors + 1;
        end
        if (slot(20) !== 25) begin
            $error("[FIR window] FIR[0] high slot[20]=%0d, expected 25", slot(20));
            errors = errors + 1;
        end
        $display("[FIR windows] FIR[0]: slot[0..20]=%0d..%0d  FIR[7]: slot[7..27]=%0d..%0d",
                 slot(0), slot(20), slot(7), slot(27));

        // ────────────────────────────────────────────────────────────────────
        // 8. MID-STREAM RESET
        // ────────────────────────────────────────────────────────────────────
        rst_n = 0;
        @(posedge clk); #1;
        rst_n = 1;
        for (k = 0; k < MEM_LEN; k = k + 1) begin
            if (slot(k) !== 0) begin
                $error("[MID RESET] slot[%0d] = %0d, expected 0", k, slot(k));
                errors = errors + 1;
            end
        end
        $display("[MID RESET] done.");

        // ────────────────────────────────────────────────────────────────────
        // RESULT
        // ────────────────────────────────────────────────────────────────────
        if (errors == 0)
            $display("PASS – all checks OK.");
        else
            $display("FAIL – %0d error(s) detected.", errors);

        $finish;
    end

    // ── Waveform dump ────────────────────────────────────────────────────────
    initial begin
        $dumpfile("shift_reg_parallel_tb.vcd");
        $dumpvars(0, shift_reg_parallel_tb);
    end

endmodule
