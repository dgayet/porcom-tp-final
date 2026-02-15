module adaptation_engine #(
    parameter NB_I      = 18,
    parameter NBF_I     = 15,
    parameter FFE_LEN   = 21,
    parameter NB        = 8,
    parameter NBF       = 7,
    parameter NBF_COEFF = 23,    // Fractional bits of cma_r input
    parameter NB_MU     = 16
)
(
    input clk,
    input rst_n,
    input enable,
    
    // Input signals from signal path
    input [NB_I-1:0]                    i_fir_out,    // filtered symbol
    input signed [NB_I*FFE_LEN-1:0]     mem_in_data,  // concatenated x[n], x[n-1], ..., x[n-FFE_LEN+1]

    // Adaption Config
    input signed [NB-1:0]               cma_r,  // CMA radius (e.g., 0.5125)
    input signed [NB_MU-1:0]            mu_cma,
    input signed [NB_MU-1:0]            mu_lms,
    
    // Weight memory interface
    input signed [FFE_LEN*NB-1:0]       coeff_flat,   // Current coefficients read from memory
    output reg signed [FFE_LEN*NB-1:0]      o_new_coeff,
    output reg signed                       o_update_en,
    
    // Status outputs
    output [31:0] iteration_count,
    output [2:0] adaptation_phase   // 0: startup, 1: CMA, 2: LMS
);

    // Config Parameters
    localparam STARTUP_DELAY = 3*FFE_LEN;  // Cycles to wait before starting adaptation. min: FFE_LEN for CMA to fill, plus some margin.
    localparam CMA_DURATION = 1000;            // Number of cycles to run CMA before switching to LMS (example value, adjust as needed)

    // Internal signals
    wire cma_enable, lms_enable;
    wire [FFE_LEN*NB-1:0] cma_weight_out;
    wire [FFE_LEN*NB-1:0] lms_weight_out;
    
    wire [2:0] phase_reg;
    
    // Instantiate adaptation controller (FSM)
    adaptation_controller adapt_ctrl (
        .clk(clk),
        .rst_n(rst_n),
        .enable(enable),
        .startup_delay(STARTUP_DELAY),
        .cma_duration(CMA_DURATION),
        .iteration_count(iteration_count),
        .adaptation_phase(phase_reg)
    );
    
    // Instantiate CMA adapter (always calculating)
    cma_full  #(
        .NB_I(NB_I),
        .NBF_I(NBF_I),
        .NB(NB),
        .NBF(NBF),
        .NBF_COEFF(NBF_COEFF),
        .NB_MU(NB_MU)
    )cma_inst(
        .i_clock(clk),
        .i_reset(rst_n),
        .i_fir_out(i_fir_out),
        .cma_r(cma_r),
        .i_xk_flat(mem_in_data),
        .i_coeff_flat(coeff_flat),
        .i_mu(mu_cma),
        .o_new_coeff(cma_weight_out)
    );
    
    // Instantiate LMS adapter (always calculating)
    // lms_adapter lms_inst (
    //     .clk(clk),
    //     .i_reset(rst_n),
    //     .i_fir_out(i_fir_out),
    //     .i_coeff_flat(coeff_flat),
    //     .i_mu(mu_lms),
    //     .weight_data_out(lms_weight_out),
    //     .weight_addr(lms_weight_addr)
    // );
    
    // Weight commit logic: route CMA or LMS output to memory
    always @(*) begin
        case (phase_reg)
            3'b000: begin  // Startup phase - no writes
                o_new_coeff = coeff_flat;  // Hold current coeffs
                o_update_en = 1'b0;
            end
            3'b001: begin  // CMA phase - commit CMA weights
                o_new_coeff = cma_weight_out;
                o_update_en = 1'b1;
            end
            3'b010: begin  // LMS phase - commit LMS weights
            // for now, use cma weights as placeholder since LMS not implemented
                o_new_coeff = cma_weight_out; // Replace with lms_weight_out when LMS is implemented
                //o_new_coeff <= lms_weight_out;
                o_update_en = 1'b1;
            end
            default: begin
                o_update_en = 1'b0;
                o_new_coeff = coeff_flat;  // Hold current coeffs
            end
        endcase
    end
    
    assign adaptation_phase = phase_reg;
    
endmodule