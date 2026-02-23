#%% IMPORTS AND PARAMETERS
import numpy as np
from ffe_func import *
from scipy.special import erfc
from tool._fixedInt import *
import numpy.random as rng
import numpy.random as rng

# Q-function
def Qfun(x):
    return 0.5 * erfc(x / np.sqrt(2))

## FIXED POINT HELPER FUNCTIONS
# Quantization function
def Q(x, W, F, signed='S', rnd='trunc', sat='saturate'):
    y = DeFixedInt(W, F, signedMode=signed, roundMode=rnd, saturateMode=sat)
    if isinstance(x, DeFixedInt):
        y.value = x.fValue
    else:
        y.value = float(x)
    return y

def FIR_fx(samples_fx, coeffs_fx, ACC_W=52, ACC_F=38):
    acc = DeFixedInt(ACC_W, ACC_F, roundMode='trunc', saturateMode='saturate')
    acc.value = 0.0
    for k in range(len(samples_fx)):
        acc = Q(acc + (samples_fx[k] * coeffs_fx[k]), ACC_W, ACC_F)
    return acc

def LMS_fx(ffe_fx, samples_fx, error_fx, mu_fx, W_W=28, W_F=23):
    for k in range(len(ffe_fx)):
        upd = samples_fx[k] * error_fx * mu_fx
        ffe_fx[k] = Q(ffe_fx[k] - upd, W_W, W_F)
    return ffe_fx

def CMA(ffe_fx, samples_fx, yk_fx, mu_fx, W_W=28, W_F=23):
    R = DeFixedInt(W_W, W_F, roundMode='trunc', saturateMode='saturate')
    a = np.array([-0.75, -0.25, 0.25, 0.75])
    # R.value = 1.64
    R.value = np.mean(np.abs(a)**4) / np.mean(np.abs(a)**2)
    error_fx = Q(yk_fx*yk_fx - R, W_W, W_F)
    for k in range(len(ffe_fx)):
        upd = samples_fx[k] * error_fx * mu_fx * yk_fx
        ffe_fx[k] = Q(ffe_fx[k] - upd, W_W, W_F)
    return ffe_fx

# thr1 = 2/norm --> punto medio entre 1 y 3
# lvl1 = 1/norm
# lvl3 = 3/norm
def slicer_fx(x_fx, thr1=0.894427191, lvl1=0.4472135955, lvl3=1.3416407865):
    t0 = Q(0.0, x_fx.width, x_fx.fractWidth)
    t1 = Q(thr1, x_fx.width, x_fx.fractWidth)
    m1 = Q(-thr1, x_fx.width, x_fx.fractWidth)

    if x_fx < m1:
        return Q(-lvl3, x_fx.width, x_fx.fractWidth)
    elif x_fx < t0:
        return Q(-lvl1, x_fx.width, x_fx.fractWidth)
    elif x_fx < t1:
        return Q(+lvl1, x_fx.width, x_fx.fractWidth)
    else:
        return Q(+lvl3, x_fx.width, x_fx.fractWidth)

def fValue_slicer(x : DeFixedInt, norm):
    x = x.fValue

    if (x>1.341):
        return Q(3/norm, 18, 15).fValue
    elif (x>0.447):
        return Q(1/norm, 18, 15).fValue
    elif (x>-1.341):
        return Q(-1/norm, 18, 15).fValue
    else:
        return Q(-3/norm, 18, 15).fValue

# BER calculation   
gray_bits = {
    -3: (0,0),
    -1: (0,1),
     1: (1,1),
     3: (1,0),
}

def sym_to_bits(sym_norm, norm):
    lvl = int(np.round(sym_norm * norm))  # should be -3,-1,1,3
    return gray_bits[lvl]

def compute_parallel_fir_outputs(mem_in_data, coeffs_fx, paralelism, ffe_len, ACC_W, ACC_F, X_W, X_F):
    """
    Compute FIR outputs for PARALELISM samples in parallel.
    mem_in_data contains FFE_LEN + PARALELISM - 1 samples.
    Returns an array of PARALELISM FIR outputs in the same order as the input samples.
    
    Memory layout (mem_in_data[0] is newest):
    - mem_in_data[0] = newest sample (latest in current batch)
    - mem_in_data[paralelism-1] = oldest sample in current batch
    - mem_in_data[paralelism:] = older samples from previous cycles
    
    Output order matches input sample order:
    - fir_outs[0] = output for oldest sample in batch
    - fir_outs[paralelism-1] = output for newest sample in batch
    """
    fir_outs = []
    for i in range(paralelism):
        # Extract the FFE_LEN samples for this output
        start_idx = paralelism - 1 - i
        samples_for_fir = mem_in_data[start_idx : start_idx + ffe_len]
        # Compute FIR
        out_fir = FIR_fx(samples_for_fir, coeffs_fx, ACC_W, ACC_F)
        out_fir = Q(out_fir, X_W, X_F)
        fir_outs.append(out_fir)
    return fir_outs

#%%
# MODULATION
PAM = 4
SR = 4e9          # symbol rate (baud)
BR = SR * np.log2(PAM)   # bit rate

BW = SR/2 # Nyquist BandWidth

print(f"- - - - - - - - - - - - - - - - - - - - -")
print(f"STARTING FP SIMULATION: PAM {PAM} EQUALIZER")
print(f"\tBit Rate: {BR/1e9:.2f}Gbps\n\tSymbol Rate: {SR/1e9:.2f}GBd\n\tBand Width: {BW/1e9:.2f}GHz\n")

# PAM-4 average symbol energy = (M^2 - 1) / 3 = 5
Es_pam = (PAM**2 - 1) / 3
norm = np.sqrt(Es_pam)   # normalizado a Es = 1

# FIXED POINT PARAMETERS
X_W, X_F = 18, 15
W_W, W_F = 28, 23
ACC_W, ACC_F = 52, 38
MU_W, MU_F = 16, 15

# SYMBOL GENERATION vs READING FROM MEM
symbol_gen = 'read' # 'gen' or 'read'

#%% SYMBOL GENERATION 

if symbol_gen == 'gen': # random symbol genearation
    print('Generating PAM4 symbols..')   
    n_symbols = 600_000*8
    PAM = 4

    # Generate PAM4 symbols: {-3, -1, +1, +3}
    symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
    # Normalize to {-0.75, -0.25, +0.25, +0.75}
    symbols = symbols_raw * 0.25

    # APPLY CHANNEl
    b = channel_fir_nyquist_loss( nyq_loss_db=114, NTAPS=11, plt_en=True)
    a = 1

    channel_symbols = np.convolve(symbols, b, mode="full")
    channel_symbols = channel_symbols[:len(symbols)]
    n_samples = len(channel_symbols)

    sample_symb_fx = arrayFixedInt(X_W, X_F, channel_symbols)
    sample_symb_fx_arr = np.array([s.fValue for s in sample_symb_fx])

    symb_fx = arrayFixedInt(X_W, X_F, symbols)
    symb_fx_arr = np.array([s.fValue for s in symb_fx])

    x = channel_symbols.copy()
    x = arrayFixedInt(X_W, X_F, x)

    bin_symb_ch = [str(s.value) for s in sample_symb_fx]
    with open("./output_mem/channel_symbols_int.mem", "w") as f:
        for s in bin_symb_ch:
            f.write(s + "\n")

    bin_symb_ch = [s.bit() for s in sample_symb_fx]
    with open("./output_mem/channel_symbols.mem", "w") as f:
        for s in bin_symb_ch:
            f.write(s + "\n")
else:  # reading fixed point symbols (in int) from mem
    print('Reading channel symbols from mem file...')    
    channel_symbols = np.loadtxt("./output_mem/channel_symbols_int.mem", dtype=int)
    x = np.zeros(len(channel_symbols)).tolist()
    x = arrayFixedInt(X_W, X_F, x)
    for i,s in enumerate(x):
        s._setValue(int(channel_symbols[i]))
    n_samples = len(x)

#%%
# parameters
FFE_LEN = 21
PARALELISM = 8
MEM_LEN = FFE_LEN + PARALELISM - 1  # Total memory needed to hold all samples for parallel processing
CENTRAL_TAP = FFE_LEN//2
FFE = np.zeros(FFE_LEN)
FFE[CENTRAL_TAP] = 8388607/2**(W_F)

mu_ffe = 1/2**10
mu_cma = 1/2**10

CMA_COUNT = 520_000*8-8
STARTUP_DELAY = 8*8-8

# quantization of variables
mu_fx = Q(mu_ffe, MU_W, MU_F)
mu_cma = Q(mu_cma, MU_W, MU_F)
FFE_fx = arrayFixedInt(W_W, W_F, FFE)
mem_in_data = np.zeros(MEM_LEN).tolist()
mem_in_data = arrayFixedInt(X_W, X_F, mem_in_data)

# scopes
mem_in_data_list = []; ydec = []
out_ffe_scope = []; cma_error_scope = []
FFE_history = []
error_scope = []

# CMA DEBUGGING
# ffe_updates = []; updates = []
# R = DeFixedInt(W_W, W_F, roundMode='trunc', saturateMode='saturate')
# aux = np.array([-0.75, -0.25, 0.25, 0.75])
# # R.value = 1.64
# R.value = np.mean(np.abs(aux)**4) / np.mean(np.abs(aux)**2)


for cycle, i in enumerate(range(0, len(x), PARALELISM)):
    if (cycle % 10000 == 0):
        print(f"cycle {cycle} (sample index {i}/{n_samples})")
    
    # Extract PARALELISM samples for this cycle
    samples_batch = x[i:i+PARALELISM]

    #print(f"Batch samples (cycle {cycle}): {[s.value for s in samples_batch]}")  # Debug: print batch samples
    # Shift memory and add new samples
    # Move existing samples to the right to make room for new samples at the front
    for k in range(MEM_LEN - 1, PARALELISM - 1, -1):
        mem_in_data[k].assign(mem_in_data[k-PARALELISM])
    
    # Insert new samples at the front in reverse order (newest at index 0)
    # samples_batch = [oldest, ..., newest] so we reverse it
    for p in range(PARALELISM):
        mem_in_data[p].assign(samples_batch[PARALELISM - 1 - p])
    
    # Store memory state
    mem_in_data_list.append([s.value for s in mem_in_data])
    
    if i == STARTUP_DELAY:
        print(f"At startup delay (cycle {cycle}) (iteration [{i}:{i+PARALELISM})):")
        print(mem_in_data)
        mem_in_data_first = [s.fValue for s in mem_in_data]
    
    # Compute FIR outputs for all PARALELISM samples in parallel
    fir_outs = compute_parallel_fir_outputs(mem_in_data, FFE_fx, PARALELISM, FFE_LEN, ACC_W, ACC_F, X_W, X_F)
    
    # Apply slicer to all outputs
    slicer_outs = [slicer_fx(out, thr1=0.50, lvl1=0.25, lvl3=0.75) for out in fir_outs]
    
    # Store all FIR outputs and slicer outputs
    for fir_out, slc_out in zip(fir_outs, slicer_outs):
        out_ffe_scope.append(fir_out)
        ydec.append(slc_out)
    
    # Weight update: only compute error for the NEWEST (last) sample
    out_ffe_latest = fir_outs[-1]
    slicer_out_latest = slicer_outs[-1]
    error_latest = out_ffe_latest - slicer_out_latest
    error_scope.append(error_latest)
    
    # Use the newest sample set [0 : FFE_LEN] for coefficient update
    samples_for_update = mem_in_data[0:FFE_LEN]
    
    if i >= STARTUP_DELAY:
        if i < CMA_COUNT + STARTUP_DELAY:
            FFE = CMA(FFE_fx, samples_for_update, out_ffe_latest, mu_cma, W_W, W_F)
        else:
            FFE = LMS_fx(FFE_fx, samples_for_update, error_latest, mu_fx, W_W, W_F)
    
    FFE_history.append(np.array([s for s in FFE_fx]))

# %% READ FROM VERILOG SIMULATION OUTPUT
rtl_signal = np.loadtxt("C:\\Users\\denis\\Documents\\beca\\porcom-tp-final\\verilog\\testbench\\out_ffe.txt", dtype=int)
rtl_signal = rtl_signal[1:]

golden_signal = np.array([s.value for s in out_ffe_scope]) # change signal according to input file

if rtl_signal.shape[0] < golden_signal.shape[0]:
    golden_signal = golden_signal[0:len(rtl_signal)]
else:
    rtl_signal = rtl_signal[0:len(golden_signal)]

print("Max abs diff:", np.max(np.abs(rtl_signal - golden_signal)))
if np.max(np.abs(rtl_signal - golden_signal)) == 0:
    print("All values match!")
else:
    for i in range(len(rtl_signal)):
        if rtl_signal[i] != golden_signal[i]:
            print(f"Mismatch at index {i}: RTL={rtl_signal[i]}, Golden={golden_signal[i]}")

#%% READ FROM VERILOG OUTPUT: FFE TAPS HISTORY
rtl_signal = np.loadtxt("C:\\Users\\denis\\Documents\\beca\\porcom-tp-final\\verilog\\testbench\\coeff_ffe.txt", dtype=int)
rtl_signal = np.array([np.flip(arr) for arr in rtl_signal])
rtl_signal = rtl_signal[2:]

golden_signal = np.array([[s.value for s in f] for f in FFE_history]) # change signal according to input file

if rtl_signal.shape[0] < golden_signal.shape[0]:
    golden_signal = golden_signal[0:len(rtl_signal)]
else:
    rtl_signal = rtl_signal[0:len(golden_signal)]

print("Max abs diff:", np.max(np.abs(rtl_signal - golden_signal)))
if np.max(np.abs(rtl_signal - golden_signal)) == 0:
    print("All values match!")
else:
    for i in range(len(rtl_signal)):
        if rtl_signal[i] != golden_signal[i]:
            print(f"Mismatch at index {i}: RTL={rtl_signal[i]}, Golden={golden_signal[i]}")

# %%
