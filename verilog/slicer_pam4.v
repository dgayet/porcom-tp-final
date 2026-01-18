module slicer_pam4
#(
    parameter NB    = 8, // word width
    parameter NF    = 6  // fraction width
)
(
    input                   i_clock,          // Clock signal
    input                   i_reset,          // Asynchronous reset
    input                   i_enable,         // Shift enable
    input                   i_valid,          // shift valid
    input  signed [NB-1:0]  i_sample,           // input sample 
    output signed [NB-1:0]  o_slicer,          // slicer output
    output signed [1:0]     o_gray_level        // gray code level
);

// Slicer for PAM4 
// normalized levels: -3/sqrt(5), -1/sqrt(5), 1/sqrt(5), 3/sqrt(5)
// decision thresholds at -2/sqrt(5), 0, 2/sqrt(5)
// fixed-point scale (S(NB, NB-1))
localparam integer SCALE = 1 << NF;

// PAM4 levels
localparam signed [NB-1:0] PAM4_P3 =  ( 3 * SCALE ) / $sqrt(5.0);
localparam signed [NB-1:0] PAM4_P1 =  ( 1 * SCALE ) / $sqrt(5.0);
localparam signed [NB-1:0] PAM4_M1 = -( 1 * SCALE ) / $sqrt(5.0);
localparam signed [NB-1:0] PAM4_M3 = -( 3 * SCALE ) / $sqrt(5.0);

// Decision thresholds
localparam signed [NB-1:0] TH_P2 =  ( 2 * SCALE ) / $sqrt(5.0);
localparam signed [NB-1:0] TH_M2 = -( 2 * SCALE ) / $sqrt(5.0);

wire signed [NB-1:0]    data_out;
wire        [1:0]       slicer_level;
always @(posedege i_clock or negedge i_reset) begin
    if (!i_reset) begin
        data_out <= {NB{1'b0}};
    end else if (i_enable & i_valid) begin
        if (i_sample >= TH_P2) begin
            data_out <= PAM4_P3;
            slicer_level <= 2'b11;
        end else if (i_sample >= 0) begin
            data_out <= PAM4_P1;
            slicer_level <= 2'b10;
        end else if (i_sample >= TH_M2) begin
            data_out <= PAM4_M1;
            slicer_level <= 2'b01;
        end else begin
            data_out <= PAM4_M3;
            slicer_level <= 2'b00;
        end
    end
end;

assign o_slicer = data_out;
assign o_gray_level = slicer_level;

endmodule;