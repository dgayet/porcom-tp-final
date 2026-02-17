module slicer_pam4
#(
    parameter NB = 8,
    parameter NBF = 7
)
(
    input                   i_clock,
    input                   i_reset,
    input                   i_enable,
    input                   i_valid,
    input   signed [NB-1:0] i_sample,
    output  signed [NB-1:0] o_slicer,
    output  [1:0]           o_gray_level
);

// Fixed-point thresholds
localparam signed [NB-1:0] TH_NEG_HALF = - (1 <<< (NBF-1)); // -0.5
localparam signed [NB-1:0] TH_POS_HALF =   (1 <<< (NBF-1)); // +0.5

// Output levels (fixed-point)
localparam signed [NB-1:0] LVL_NEG_075 = - (3 <<< (NBF-2)); // -0.75
localparam signed [NB-1:0] LVL_NEG_025 = - (1 <<< (NBF-2)); // -0.25
localparam signed [NB-1:0] LVL_POS_025 =   (1 <<< (NBF-2)); // +0.25
localparam signed [NB-1:0] LVL_POS_075 =   (3 <<< (NBF-2)); // +0.75

// Combinational logic: no delay
assign o_slicer = (i_sample < TH_NEG_HALF) ? LVL_NEG_075 :
                  (i_sample < 0)            ? LVL_NEG_025 :
                  (i_sample < TH_POS_HALF)  ? LVL_POS_025 :
                                              LVL_POS_075;

assign o_gray_level = (i_sample < TH_NEG_HALF) ? 2'b00 :
                      (i_sample < 0)            ? 2'b01 :
                      (i_sample < TH_POS_HALF)  ? 2'b11 :
                                                  2'b10;

endmodule