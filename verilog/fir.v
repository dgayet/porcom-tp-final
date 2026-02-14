module fir
    #(
        parameter FIR_LEN   = 21,
        parameter NB_COEFF  = 28,
        parameter NBF_COEFF = 23,
        parameter NB_IN     = 18,
        parameter NBF_IN    = 15,
        parameter NB_OUT    = 18,
        parameter NBF_OUT   = 15
        )
    (
        output signed [NB_OUT-1:0]              o_sample, // NB: NBF + 1 (sign bit)
        input                                   clk,
        input                                   i_reset,
        input  signed [FIR_LEN*NB_IN -1:0]      i_data_reg,
        input  signed [FIR_LEN*NB_COEFF -1:0]   i_coeff,
        input                                   i_en,
        input                                   i_valid                    
    );

        localparam NB_ADD     = NB_COEFF    + NB_IN + $clog2(FIR_LEN) + 1;
        localparam NBF_ADD    = NBF_COEFF   + NBF_IN; //check
        localparam NBI_ADD    = NB_ADD      - NBF_ADD;
        localparam NBI_OUTPUT = NB_OUT      - NBF_OUT;
        localparam NB_SAT     = (NBI_ADD)   - (NBI_OUTPUT);

        
        reg  signed [NB_IN + NB_COEFF -1 : 0]   prod        [FIR_LEN-1 : 0]; //! Partial Products
        reg  signed [NB_ADD - 1 : 0]            sum;

        // unpacking
        reg signed [NB_IN-1:0]      xk    [0:FIR_LEN-1];
        reg signed [NB_COEFF-1:0]   coeff [0:FIR_LEN-1];

        integer k;
        always @(*) begin
            for (k = 0; k < FIR_LEN; k = k + 1) begin
                xk[k] = i_data_reg[k*NB_IN +: NB_IN];
                coeff[k]  = i_coeff[k*NB_COEFF +: NB_COEFF];
            end
        end


        integer ptr;
        always @(*) begin
            if (!i_reset) begin
                for (ptr = 0; ptr < FIR_LEN; ptr = ptr + 1)
                    prod[ptr] =  {(NB_IN + NB_COEFF){1'b0}};
            end
            else if (i_en && i_valid) begin
                for(ptr=0;ptr<FIR_LEN;ptr=ptr+1) begin:mult
                        prod[ptr] = coeff[ptr] * xk[ptr];
                    end  
          end  
        end
    
        integer ptr3;
        always @(*) begin:accum
        sum = {NB_ADD{1'b0}};
        for(ptr3=0;ptr3<FIR_LEN;ptr3=ptr3+1) begin:adder 
            sum = sum + prod[ptr3];
        end
        end
        
        // Output
        wire of_pos =  ~|sum[NB_ADD-1 -: NB_SAT+1];
        wire of_neg =  &sum[NB_ADD-1 -: NB_SAT+1];
        assign o_sample = ( of_pos || of_neg) ? sum[NB_ADD-(NBI_ADD-NBI_OUTPUT) - 1 -: NB_OUT] :
                    (sum[NB_ADD-1]) ? {{1'b1},{NB_OUT-1{1'b0}}} : {{1'b0},{NB_OUT-1{1'b1}}};

endmodule