// Convention:
// h[N-1] = newest tap, multiplies x[N-1] (current symbol)
// h[k] multiplies x[k] where index 0 is oldest sample
// Matches Python CMA implementation
//
// Parallelism:
// PARALLELISM = 8 samples per clock
// MEM_LEN = FIR_LEN + PARALLELISM - 1 = 28 samples stored
// Input samples ordered: i_sample[NB_IN-1 -: NB_IN] = oldest, i_sample[PARALLELISM*NB_IN-1 -: NB_IN] = newest
// FIR[i] window offset i (i=0 oldest output, i=7 newest output)
// Adaptation engine fed only newest sample's context (from FIR[7] and mem_flat[MEM_LEN*NB_IN-1 -: FIR_LEN*NB_IN])

module top
#(
    parameter FIR_LEN   = 21,
    parameter NB_COEFF  = 28,
    parameter NBF_COEFF = 23,
    parameter NB_IN     = 18,
    parameter NBF_IN    = 15,
    parameter NB_OUT    = 18,
    parameter NBF_OUT   = 15,
    parameter NB_MU     = 16,
    parameter PARALLELISM = 8
)
(
    input                                   i_clock,
    input                                   i_reset,
    input                                   i_en,
    input                                   i_valid,
    input  signed [PARALLELISM*NB_IN-1:0]   i_sample,
    input  signed [NB_MU-1:0]               i_mu,
    output signed [PARALLELISM*NB_OUT-1:0]  o_sample
);
    localparam MEM_LEN = FIR_LEN + PARALLELISM - 1;  // 28
    localparam signed [NB_COEFF-1:0] CMA_R = 28'sd4299161; // 0.5125

    wire signed [FIR_LEN*NB_COEFF-1:0]      coeff_flat;
    wire signed [FIR_LEN*NB_COEFF-1:0]      new_coeff_flat;
    wire signed [MEM_LEN*NB_IN-1:0]         mem_flat;

    wire update_en;
    wire signed [PARALLELISM*NB_OUT-1:0]    fir_out_flat;
    wire signed [PARALLELISM*NB_OUT-1:0]    slicer_out_flat;

    // ── SHIFT REGISTER: stores MEM_LEN = 28 samples, accepts PARALLELISM = 8 per cycle ──
    shift_register #(
        .N          (FIR_LEN),
        .NB         (NB_IN),
        .PARALLELISM(PARALLELISM)
    ) u_sample_line (
        .i_clock (i_clock),
        .i_reset (i_reset),
        .i_enable(i_en),
        .i_valid (i_valid),
        .i_data  (i_sample),
        .o_data  (mem_flat)
    );

    // ── 8 PARALLEL FIR INSTANCES ──────────────────────────────────────────────────────
    // FIR[i] processes samples at offset i (i=0 oldest output, i=7 newest output)
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
                .clk       (i_clock),
                .i_reset   (i_reset),
                .i_data_reg(mem_flat[(fir_i+FIR_LEN)*NB_IN - 1 -: FIR_LEN*NB_IN]),
                .i_coeff   (coeff_flat),
                .i_en      (i_en),
                .i_valid   (i_valid)
            );
        end
    endgenerate

    // ── 8 PARALLEL SLICER INSTANCES ───────────────────────────────────────────────────
    genvar slicer_i;
    generate
        for (slicer_i = 0; slicer_i < PARALLELISM; slicer_i = slicer_i + 1) begin : slicer_gen
            slicer_pam4 #(
                .NB (NB_OUT),
                .NBF(NBF_OUT)
            ) slicer_inst (
                .i_clock   (i_clock),
                .i_reset   (i_reset),
                .i_enable  (i_en),
                .i_valid   (i_valid),
                .i_sample  (fir_out_flat[(slicer_i+1)*NB_OUT - 1 -: NB_OUT]),
                .o_slicer  (slicer_out_flat[(slicer_i+1)*NB_OUT - 1 -: NB_OUT]),
                .o_gray_level()
            );
        end
    endgenerate

    // ── ADAPTATION ENGINE: fed newest sample's context only ────────────────────────────
    // Note: FIR is combinational, so shift_reg output and FIR output are already aligned
    adaptation_engine #(
        .NB_I      (NB_OUT),
        .NBF_I     (NBF_OUT),
        .FFE_LEN   (FIR_LEN),
        .NB        (NB_COEFF),
        .NBF       (NBF_COEFF),
        .NBF_COEFF (NBF_COEFF),
        .NB_MU     (NB_MU)
    ) DUT (
        .clk         (i_clock),
        .rst_n       (i_reset),
        .enable      (i_en),
        .i_fir_out   (fir_out_flat[PARALLELISM*NB_OUT - 1 -: NB_OUT]),       // newest FIR output
        .i_slicer_out(slicer_out_flat[PARALLELISM*NB_OUT - 1 -: NB_OUT]),   // newest slicer output
        .mem_in_data (mem_flat[MEM_LEN*NB_IN - 1 -: FIR_LEN*NB_IN]),        // newest 21-sample window
        .cma_r       (CMA_R),
        .mu_cma      (i_mu),
        .mu_lms      (i_mu),
        .coeff_flat  (coeff_flat),
        .o_new_coeff (new_coeff_flat),
        .o_update_en (update_en)
    );
    
    // ── WEIGHTS MEMORY ────────────────────────────────────────────────────────────────
    weights_mem #(
        .FIR_LEN    (FIR_LEN),
        .NB_COEFF   (NB_COEFF),
        .NBF_COEFF  (NBF_COEFF),
        .CENTRAL_TAP(10)
    ) FFE_MEM (
        .i_clock     (i_clock),
        .i_reset     (i_reset),
        .i_update_en (update_en),
        .i_w_new_flat(new_coeff_flat),
        .o_w_flat    (coeff_flat)
    );

    // ── OUTPUT: pack all 8 FIR outputs ────────────────────────────────────────────────
    assign o_sample = fir_out_flat;

endmodule