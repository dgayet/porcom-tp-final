module adaptation_controller (
    input clk,
    input rst_n,
    input enable,
    
    // Configuration
    input [31:0] startup_delay,
    input [31:0] cma_duration,
    
    // Outputs
    output reg [31:0] iteration_count,
    output reg [2:0] adaptation_phase  // 0: startup, 1: CMA, 2: LMS
);

    // FSM states
    localparam STARTUP = 3'b000;
    localparam CMA = 3'b001;
    localparam LMS = 3'b010;
    
    reg [31:0] counter;
    reg [2:0] state, next_state;
    
    // FSM state transitions
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= STARTUP;
        end
        else begin
            state <= next_state;
        end
    end
    
    always @(*) begin
        next_state = state;
        case (state)
            STARTUP: begin
                if (counter >= startup_delay) begin
                    next_state = CMA;
                end
            end
            CMA: begin
                if (counter >= (startup_delay + cma_duration)) begin
                    next_state = LMS;
                end
            end
            LMS: begin
                next_state = LMS;  // Stay in LMS indefinitely
            end
        endcase
    end
    
    // output logic
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            iteration_count <= 32'b0;
            adaptation_phase <= STARTUP;
        end
        else if (enable) begin
            adaptation_phase <= next_state;
            
            // Calculate iteration count based on phase
            case (state)
                STARTUP: begin
                    iteration_count <= 32'b0;
                end
                CMA: begin
                    iteration_count <= counter - startup_delay;
                end
                LMS: begin
                    iteration_count <= counter - startup_delay - cma_duration;
                end
            endcase
        end
    end

    // counter logic
    always @(posedge clk or negedge rst_n) begin // add valid and enable logic ..
        if (!rst_n) begin
            counter <= 32'b0;
        end else if (enable) begin
            counter <= counter + 1'b1;
        end
    end
    
endmodule