module cma_update
#(
    parameter                   NB_I = 18,
    parameter                   NBF_I = 15,
    parameter                   FFE_LEN =21,
    parameter                   NB = 8,
    parameter                   NBF = 7,
    parameter                   NB_MU = 16
)
(
    input                                   i_clock,
    input                                   i_reset,
    input signed    [NB_I-1:0]              i_fir_out,
    input signed    [NB-1:0]                cma_r,  // 0.5125 
    //input signed    [NB_I*(FFE_LEN-1)-1:0]  i_shift_reg,
    //input signed    [NB_MU-1:0]             i_mu,
    output signed   [NB-1:0]                cma_error
);

wire signed [2*NB_I-1:0] y_sq_full;
assign y_sq_full = i_fir_out * i_fir_out;

wire signed [2*NB_I-1:0] R_full;
assign R_full = cma_r <<< (2*NBF_I - NBF);

// error = y^2 - R
wire signed [2*NB_I:0] err_full;
assign err_full = y_sq_full - R_full;

assign cma_error = err_full >>> (2*NBF_I - NBF);

endmodule;