module delay_line
#(
    parameter WIDTH = 168,  // Data width
    parameter DEPTH = 2     // Number of pipeline stages
)
(
    input                       i_clock,
    input                       i_reset,
    input                       i_valid,
    input  signed [WIDTH-1:0]   i_data,
    output signed [WIDTH-1:0]   o_data
);

    reg signed [WIDTH-1:0] pipeline [0:DEPTH-1];
    integer i;

    always @(posedge i_clock or negedge i_reset) begin
        if (!i_reset) begin
            for (i = 0; i < DEPTH; i = i + 1)
                pipeline[i] <= {WIDTH{1'b0}};
        end
        else if (i_valid) begin
            pipeline[0] <= i_data;
            for (i = 1; i < DEPTH; i = i + 1)
                pipeline[i] <= pipeline[i-1];
        end
    end

    assign o_data = pipeline[DEPTH-1];

endmodule