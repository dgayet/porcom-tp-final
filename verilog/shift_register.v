module shift_register 
#(
    parameter N = 21, // length of shift register
    parameter NB = 18 // word width
)
(
    input               i_clock,          // Clock signal
    input               i_reset,          // Asynchronous reset
    input               i_enable,         // Shift enable
    input               i_valid,          // shift valid
    input  [NB-1:0]     i_data,           // new bit 
    output [N*NB-1:0]   o_data          // Shift register output
);


reg [N*NB-1:0]  shift_reg;
always @(posedge i_clock or negedge i_reset) begin
    if (!i_reset) begin
        shift_reg <= {N*NB{1'b0}}; // Clear register on reset
    end else if (i_enable & i_valid) begin
        // shift RIGHT, insert newest at MSB
        shift_reg <= {i_data, shift_reg[N*NB-1:NB]};
    end
end

assign o_data = shift_reg;

endmodule
