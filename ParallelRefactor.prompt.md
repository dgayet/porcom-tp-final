# FFE Parallel Refactor Plan (PARALLELISM=8)

## Overview

Process 8 input samples per clock by:
- Extending `shift_register.v` with a `PARALLELISM` parameter
- Refactoring `top.v` to instantiate 8 parallel FIR + slicer pairs
- Feeding the adaptation engine exclusively the newest sample's FIR output and its corresponding 21-sample window (no changes to adaptation engine internals)
- Correcting FSM timing localparams in `adaptation_engine.v` to match the Python simulation cycle counts
- Creating a new `parallel_top_tb.v` testbench

---

## Step 1 — Extend `shift_register.v`

Add `PARALLELISM` parameter (default `1` for backward compatibility).

- Stored depth becomes `MEM_LEN = N + PARALLELISM - 1` (e.g. 21 + 8 - 1 = 28)
- Output width: `MEM_LEN * NB` bits
- Input bus: `i_data[PARALLELISM*NB-1:0]`
  - `i_data[0 +: NB]` = oldest sample of batch
  - `i_data[(PARALLELISM-1)*NB +: NB]` = newest sample of batch
- On each `enable & valid`, shift out `PARALLELISM` oldest samples and shift in `PARALLELISM` new ones:

```verilog
shift_reg <= { i_data[PARALLELISM*NB-1:0], shift_reg[MEM_LEN*NB-1 : PARALLELISM*NB] };
```

- Output `o_data[0 +: NB]` = globally oldest sample, `o_data[(MEM_LEN-1)*NB +: NB]` = globally newest

- create testbench `shift_register_parallel_tb.v` to verify correct shifting and output ordering with `PARALLELISM=8`
---

## Step 2 — Refactor `top.v`

### New parameters / localparams

```verilog
parameter PARALLELISM = 8
localparam MEM_LEN = FIR_LEN + PARALLELISM - 1  // = 28
```

### I/O changes

- `i_sample` → `input signed [PARALLELISM*NB_IN-1:0] i_sample`
  - `i_sample[0 +: NB_IN]` = oldest in batch, `i_sample[(PARALLELISM-1)*NB_IN +: NB_IN]` = newest
- `o_sample` → `output signed [PARALLELISM*NB_OUT-1:0] o_sample`
  - Same ordering: index 0 = oldest output, index PARALLELISM-1 = newest output

### shift_register instantiation

```verilog
shift_register #(.N(MEM_LEN), .NB(NB_IN), .PARALLELISM(PARALLELISM)) u_sample_line (
    .i_clock(i_clock), .i_reset(i_reset),
    .i_enable(i_en), .i_valid(i_valid),
    .i_data(i_sample),      // PARALLELISM*NB_IN wide
    .o_data(mem_flat)       // MEM_LEN*NB_IN wide
);
```

### 8 parallel FIR instances (generate)

FIR instance `i` (i=0 oldest output, i=7 newest output):
- `i_data_reg = mem_flat[i*NB_IN +: FIR_LEN*NB_IN]`
- All share the same `coeff_flat` from `weights_mem`

```verilog
wire signed [PARALLELISM*NB_OUT-1:0] fir_out_flat;

genvar i;
generate
  for (i = 0; i < PARALLELISM; i = i + 1) begin : fir_gen
    fir #(...) fir_inst (
      .i_data_reg(mem_flat[i*NB_IN +: FIR_LEN*NB_IN]),
      .i_coeff(coeff_flat),
      .o_sample(fir_out_flat[i*NB_OUT +: NB_OUT]),
      ...
    );
  end
endgenerate
```

### 8 parallel slicer instances (generate)

```verilog
wire signed [PARALLELISM*NB_OUT-1:0] slicer_out_flat;

genvar j;
generate
  for (j = 0; j < PARALLELISM; j = j + 1) begin : slicer_gen
    slicer_pam4 #(...) slicer_inst (
      .i_sample(fir_out_flat[j*NB_OUT +: NB_OUT]),
      .o_slicer(slicer_out_flat[j*NB_OUT +: NB_OUT]),
      ...
    );
  end
endgenerate
```

### Adaptation engine connections (unchanged interface)

Feed only the **newest sample's** context (top slice of mem_flat):

```verilog
adaptation_engine #(...) DUT (
    .i_fir_out    ( fir_out_flat   [PARALLELISM*NB_OUT-1 -: NB_OUT]   ),  // newest FIR output
    .i_slicer_out ( slicer_out_flat[PARALLELISM*NB_OUT-1 -: NB_OUT]   ),  // newest slicer output
    .mem_in_data  ( mem_flat       [MEM_LEN*NB_IN-1      -: FIR_LEN*NB_IN] ),  // newest 21-sample window
    ...
);
```

Note: `-:` anchors at the MSB (newest-sample end) and steps downward — consistent across all three signals and matches the Python memory layout `mem_in_data[0:FFE_LEN]` = newest 21 samples.

### Output assignment

```verilog
assign o_sample = fir_out_flat;
```

---

## Step 3 — Fix FSM timing in `adaptation_engine.v`

**Derivation:**
- Python `STARTUP_DELAY` = 56 samples → 56 / 8 = **7 cycles**
- Python `CMA_DURATION` = 4,159,992 samples → 4,159,992 / 8 = **519,999 cycles**

```verilog
// Before
localparam STARTUP_DELAY = 3*FFE_LEN;   // = 63
localparam CMA_DURATION  = 500000;

// After
localparam STARTUP_DELAY = 7;
localparam CMA_DURATION  = 519999;
```

All other logic in `adaptation_engine.v`, `cma_adapter.v`, `lms_adapter.v`, `cma_base.v`, `lms_base.v` is **unchanged**.

---

## Step 4 — New testbench `testbench/parallel_top_tb.v`

- Same parameters as `cma_tb.v`: `FIR_LEN=21`, `NB_COEFF=28`, `NB_IN=18`, `NBF_IN=15`, `NB_MU=16`, `mu=32`
- Read all samples from `channel_symbols.mem` into a flat array
- Loop in steps of 8; pack each group of 8 into `i_sample`:
  - `i_sample[0 +: NB_IN]` = sample[n]   (oldest in batch, first in file order)
  - `i_sample[7*NB_IN +: NB_IN]` = sample[n+7] (newest in batch)
- Each posedge: capture all 8 words from `o_sample`, write as signed integers to `out_ffe_parallel.txt` (one per line, oldest→newest)
- Also dump all `FIR_LEN` coefficient values per cycle to `coeff_ffe_parallel.txt` for convergence comparison

---

## Verification

1. Run `python/fixed_point_paralel.py` → generates golden `out_ffe_parallel_golden.txt`
2. Simulate `parallel_top_tb.v` in ModelSim / Icarus Verilog
3. Diff RTL output vs Python output line-by-line (skip first `STARTUP_DELAY` cycles)
4. Check coefficient trajectory in `coeff_ffe_parallel.txt` matches Python weight history

---

## Unchanged modules

- `adaptation_engine.v` — interface unchanged, only localparams corrected
- `adaption_controller.v` — unchanged
- `cma_adapter.v` — unchanged
- `cma_base.v` — unchanged
- `lms_adapter.v` — unchanged
- `lms_base.v` — unchanged
- `slicer_pam4.v` — unchanged
- `weights_mem.v` — unchanged
- `delay_line.v` — unchanged
