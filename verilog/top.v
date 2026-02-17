// Convention:
// h[N-1] = newest tap, multiplies x[N-1] (current symbol)
// h[k] multiplies x[k] where index 0 is oldest sample
// Matches Python CMA implementation

module top
#(
    parameter FIR_LEN   = 21,
    parameter NB_COEFF  = 8,
    parameter NBF_COEFF = 7,
    parameter NB_IN     = 8,
    parameter NBF_IN    = 7,
    parameter NB_OUT    = 18,
    parameter NBF_OUT   = 17,
    parameter NB_MU     = 16
)
(
    input                       i_clock,
    input                       i_reset,
    input                       i_en,
    input                       i_valid,
    input  signed [NB_IN-1:0]   i_sample,
    input  signed [NB_MU-1:0]   i_mu,
    output signed [NB_OUT-1:0]  o_sample
);
    localparam signed [NB_COEFF-1:0] CMA_R = 28'sd4299161; // 0.5125

    wire signed [FIR_LEN*NB_COEFF-1:0]   coeff_flat;
    wire signed [FIR_LEN*NB_COEFF-1:0]   new_coeff_flat;
    wire signed [FIR_LEN*NB_IN-1:0]  xk_flat;

    wire update_en;
    wire signed [NB_OUT-1:0]  fir_out;
    wire signed [NB_OUT-1:0]  slicer_out;

    // Signal Path
    shift_register #(
        .N(FIR_LEN),
        .NB(NB_IN)
    ) u_sample_line (
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
        .o_sample(fir_out),
        .clk(i_clock),
        .i_reset(i_reset),
        .i_data_reg(xk_flat),
        .i_coeff(coeff_flat),
        .i_en(i_en),
        .i_valid(i_valid)
    );

    // Slicer
    slicer_pam4 #(
        .NB(NB_OUT),
        .NBF(NBF_OUT)
    ) slicer_inst (
        .i_clock(i_clock),
        .i_reset(i_reset),
        .i_enable(i_en),
        .i_valid(i_valid),
        .i_sample(fir_out),
        .o_slicer(slicer_out),
        .o_gray_level()
    );

    // Add pipeline delay for xk_flat to align with fir_out
    wire signed [FIR_LEN*NB_IN-1:0]  xk_flat_aligned;
    delay_line #(
        .WIDTH(FIR_LEN*NB_IN),
        .DEPTH(1)
    ) xk_delay (
        .i_clock(i_clock),
        .i_reset(i_reset),
        .i_valid(i_valid),
        .i_data(xk_flat),
        .o_data(xk_flat_aligned)
    );

    // Adaptation engine
    adaptation_engine #(
        .NB_I(NB_OUT),
        .NBF_I(NBF_OUT),
        .FFE_LEN(FIR_LEN),
        .NB(NB_COEFF),
        .NBF(NBF_COEFF),
        .NB_MU(NB_MU)
    ) DUT (
        .clk(i_clock),
        .rst_n(i_reset),
        .enable(i_en),
        .i_fir_out(fir_out),
        .i_slicer_out(slicer_out),
        .mem_in_data(xk_flat), 
        .cma_r(CMA_R),
        .mu_cma(i_mu),
        .mu_lms(i_mu),
        .coeff_flat(coeff_flat),
        .o_new_coeff(new_coeff_flat),
        .o_update_en(update_en)
        //debug
        //.cma_error(error_out)
    );
    
    // ! Weights registering
    weights_mem #(
        .FIR_LEN(FIR_LEN),
        .NB_COEFF(NB_COEFF),
        .NBF_COEFF(NBF_COEFF),
        .CENTRAL_TAP(10)
    ) FFE_MEM (
        .i_clock(i_clock),
        .i_reset(i_reset),
        .i_update_en(update_en),
        .i_w_new_flat(new_coeff_flat),
        .o_w_flat(coeff_flat)
    );




    assign o_sample = fir_out;

endmodule