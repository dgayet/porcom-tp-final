// shift_register: parallel input shift register
//
// When PARALLELISM=1 (default) the interface and behaviour are identical to
// the original serial version.
//
// When PARALLELISM=P > 1:
//   - i_data  is P*NB bits wide.
//     i_data[NB-1 : 0]         = oldest sample of the incoming batch
//     i_data[NB*P-1 -: NB]  = newest  sample of the incoming batch
//   - The stored depth grows to MEM_LEN = N + P - 1 samples -> o_data  is MEM_LEN*NB bits wide.
//     o_data[NB-1 : 0]              = globally oldest sample in memory
//     o_data[(NB*P-1 -: NB] = globally newest sample in memory

module shift_register
#(
    parameter N           = 21, // number of FIR taps (usable window depth)
    parameter NB          = 18, // word width per sample
    parameter PARALLELISM = 1   // samples accepted per clock cycle (P)
)
(
    input                                   i_clock,
    input                                   i_reset,   // active-low async reset
    input                                   i_enable,
    input                                   i_valid,
    input      [PARALLELISM*NB-1:0]         i_data,    // batch: [NB-1 : 0]=oldest, [NB*P-1 -: NB]=newest
    output     [(N+PARALLELISM-1)*NB-1:0]   o_data     // full memory snapshot
);

    localparam MEM_LEN = N + PARALLELISM - 1;

    reg [MEM_LEN*NB-1:0] shift_reg;

    always @(posedge i_clock or negedge i_reset) begin
        if (!i_reset) begin
            shift_reg <= {MEM_LEN*NB{1'b0}};
        end else if (i_enable & i_valid) begin
            // Retire the PARALLELISM oldest samples (LSBs), insert new batch at MSB end.
            // Post-shift layout (index 0 = oldest):
            //   [MEM_LEN-1] = newest arriving sample  (i_data[(P-1)*NB +: NB])
            //   [MEM_LEN-P] = oldest arriving sample  (i_data[0 +: NB])
            //   [MEM_LEN-P-1 .. 0] = survivors from previous cycle
            shift_reg <= {i_data, shift_reg[MEM_LEN*NB-1 : PARALLELISM*NB]};
        end
    end

    assign o_data = shift_reg;

endmodule
