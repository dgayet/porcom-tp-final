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

reg    [NB-1:0] data_out;
reg    [1:0]    gray_level;

// Fixed-point thresholds
localparam signed [NB-1:0] TH_NEG_HALF = - (1 <<< (NBF-1)); // -0.5
localparam signed [NB-1:0] TH_POS_HALF =   (1 <<< (NBF-1)); // +0.5

// Output levels (fixed-point)
localparam signed [NB-1:0] LVL_NEG_075 = - (3 <<< (NBF-2)); // -0.75
localparam signed [NB-1:0] LVL_NEG_025 = - (1 <<< (NBF-2)); // -0.25
localparam signed [NB-1:0] LVL_POS_025 =   (1 <<< (NBF-2)); // +0.25
localparam signed [NB-1:0] LVL_POS_075 =   (3 <<< (NBF-2)); // +0.75

always @(posedge i_clock or negedge i_reset) begin
    if (!i_reset) begin
        data_out    <= {NB{1'b0}};
        gray_level  <= {2{1'b0}};
    end else if (i_enable & i_valid) begin
        if (i_sample < TH_NEG_HALF) begin
            gray_level  <= 2'b00;
            data_out <= LVL_NEG_075;
        end
        else if (i_sample < 0) begin
            gray_level  <= 2'b01;
            data_out <= LVL_NEG_025;
        end
        else if (i_sample < TH_POS_HALF) begin
            gray_level  <= 2'b11;
            data_out <= LVL_POS_025;
        end
        else begin
            gray_level  <= 2'b10;
            data_out <= LVL_POS_075;
        end
    end
end

assign o_slicer     = data_out;
assign o_gray_level = gray_level;

endmodule