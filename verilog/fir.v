module fir
    #(
        parameter FIR_LEN   = 21,
        parameter NB_COEFF  = 8,
        parameter NBF_COEFF = 7,
        parameter NB_IN     = 8,
        parameter NBF_IN    = 7,
        parameter NBF       = 7,
        parameter NB_OUT    = 18,
        parameter NBF_OUT   = 17
        )
    (
        output signed [NB_OUT-1:0]      o_sample, // NB: NBF + 1 (sign bit)
        input                           clk,
        input                           i_reset,
        input  signed [NB_IN-1:0]       i_is_data,
        input                           i_en,
        input                           i_valid                    
    );

        localparam NB_ADD     = NB_COEFF    + NB_IN + 2;
        localparam NBF_ADD    = NBF_COEFF   + NBF_IN; //check
        localparam NBI_ADD    = NB_ADD      - NBF_ADD;
        localparam NBI_OUTPUT = NB_OUT      - NBF_OUT;
        localparam NB_SAT     = (NBI_ADD)   - (NBI_OUTPUT);


        wire signed [NB_COEFF - 1 : 0]          coeff       [FIR_LEN-1 : 0];
        reg  signed [NB_IN - 1 : 0]             register    [FIR_LEN-2 : 0]; //! Matrix for registers
        reg  signed [NB_IN + NB_COEFF -1 : 0]   prod        [FIR_LEN-1 : 0]; //! Partial Products
        reg  signed [NB_ADD - 1 : 0]            sum;


        // FFE coefficients (S(8,7), floor quantization)
        assign coeff[20] = 8'b11111111; //  -1
        assign coeff[19] = 8'b11111111; //  -1
        assign coeff[18] = 8'b11111111; //  -1
        assign coeff[17] = 8'b00000000; //   0
        assign coeff[16] = 8'b11111111; //  -1
        assign coeff[15] = 8'b00000000; //   0
        assign coeff[14] = 8'b11111111; //  -1
        assign coeff[13] = 8'b00000000; //   0
        assign coeff[12] = 8'b11111111; //  -1
        assign coeff[11] = 8'b00000000; //   0
        assign coeff[10] = 8'b01111111; //  127 (sat)
        assign coeff[9]  = 8'b10000000; // -128 (sat)
        assign coeff[8]  = 8'b11011101; //  -35
        assign coeff[7]  = 8'b11110001; //  -15
        assign coeff[6]  = 8'b11110111; //   -9
        assign coeff[5]  = 8'b11111010; //   -6
        assign coeff[4]  = 8'b00101100; //   44
        assign coeff[3]  = 8'b11110111; //   -9
        assign coeff[2]  = 8'b11111010; //   -6
        assign coeff[1]  = 8'b11111100; //   -4
        assign coeff[0]  = 8'b11111101; //   -3
        

        //! ShiftRegister model
        integer ptr1;
        integer ptr2;
        always @(posedge clk or negedge i_reset) begin:shiftRegister
            if (!i_reset) begin
            for(ptr1=0;ptr1<FIR_LEN-1;ptr1=ptr1+1) begin:init
                register[ptr1] <= {(NB_IN){1'b0}};
            end
            end else begin
            if (i_en == 1'b1) begin
                for(ptr2=0;ptr2<FIR_LEN-1;ptr2=ptr2+1) begin:srmove
                    if(ptr2==FIR_LEN-2)
                        if (i_valid)
                            register[ptr2] <= i_is_data;
                        else
                            register[ptr2] <= {(NB_IN){1'b0}};
                    else
                        register[ptr2] <= register[ptr2+1];
                    end   
                end
            end
        end

        integer ptr;
        always @(posedge clk or negedge i_reset) begin// @(*) begin
          prod[FIR_LEN-1] <= coeff[FIR_LEN-1] * i_is_data;
          for(ptr=0;ptr<FIR_LEN-1;ptr=ptr+1) begin:mult
                prod[ptr] <= coeff[ptr] * register[ptr];
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