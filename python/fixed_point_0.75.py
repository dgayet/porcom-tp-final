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
# mem_in_data: S(16,14) 
# FFE coefs: S(17,16) --> precision: 1.52587890625e-05, 
# in FFE Coefs: min number: -1e-4
# Accumulator / out_ffe / error: S(26,18) 
# mu: S(21,20) --> precision: 9.5367431640625e-07
X_W, X_F = 18, 15
W_W, W_F = 28, 23
ACC_W, ACC_F = 52, 38
MU_W, MU_F = 16, 15

#%% BASE CASE: No channel, no noise, no pulse shaping -> FFE must be identity

# symbol generation
n_symbols = 10_000
# Generate PAM4 symbols: {-3, -1, +1, +3}
symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
# Normalize to {-0.75, -0.25, +0.25, +0.75}
symbols = symbols_raw * 0.25


# Apply Channel
channel_symbols = symbols.copy()

# Noise Generation
channel_symbols = symbols.copy()

# FFE Parameters
FFE_LEN = 21
CENTRAL_TAP = FFE_LEN//2
FFE = np.zeros(FFE_LEN)
FFE[CENTRAL_TAP] = 1.0
mu_ffe = 1e-3
mu_cma = 1e-4

# quantization of variables
mu_fx = Q(mu_ffe, MU_W, MU_F)
mu_cma = Q(mu_cma, MU_W, MU_F)
sample_symb_fx = arrayFixedInt(X_W, X_F, channel_symbols)

FFE = arrayFixedInt(W_W, W_F, FFE)
mem_in_data = np.zeros(FFE_LEN).tolist()
mem_in_data = arrayFixedInt(X_W, X_F, mem_in_data)

ydec = []
FFE_history = []

ffe_out_scope    = []
in_ffe_scope     = []
error_scope      = []

STARTUP_DELAY = 3*len(FFE)
CMA_COUNT = int(1e3)
for i, sample in enumerate(sample_symb_fx):
    # Shift memory
    for k in range(FFE_LEN-1, 0, -1):
        mem_in_data[k].assign(mem_in_data[k-1])
    mem_in_data[0].assign(sample)
    in_ffe_scope.append(sample.fValue)

    # FFE output
    #out_ffe = sample
    out_ffe         = FIR_fx(mem_in_data,FFE)
    out_ffe = Q(out_ffe, X_W, X_F)
    ffe_out_scope.append(out_ffe)

    out_slicer = slicer_fx(out_ffe, thr1=0.50, lvl1=0.25, lvl3=0.75)
    # decider
    #out_slicer = slicer_fx(out_ffe)
    ydec.append(out_slicer)

    error_slicer = Q(out_ffe-out_slicer, ACC_W, ACC_F)
    error_scope.append(error_slicer.fValue)
    if i > STARTUP_DELAY:
        if i < CMA_COUNT:
            FFE = CMA(FFE,mem_in_data, out_ffe, mu_cma)
        else:
            FFE = LMS_fx(FFE,mem_in_data, error_slicer,mu_fx)
    FFE_history.append(FFE.copy())

FFE_history = [[tap.fValue for tap in ffe] for ffe in FFE_history]

#valid = np.arange(CENTRAL_TAP, n_symbols)
ydec_slic = np.array([y.fValue for y in ydec])
ydec = np.array([y.fValue for y in ydec])
in_ffe_scope = np.array(in_ffe_scope)
ffe_out_scope = np.array([s.fValue for s in ffe_out_scope])
sample_symb_fx_arr = np.array([s.fValue for s in sample_symb_fx])


print("Caso Base: el canal es la identidad, sin ruido, sin pulse shaping:" \
"espero que el FFE quede igual al impulso, solamente habrá que alinear la secuencia de entrada" \
"con la de salida del FFE (delay porque el tap central está en el medio del filtro)")

GET_SER(ydec_slic, CENTRAL_TAP, sample_symb_fx_arr, start_ber=0)
print("Central tap final =", FFE[CENTRAL_TAP])

plot_ffe(FFE_history, 10, f"FFE Coefficient Evolution - mu={mu_ffe:.1e}")


# %% AWGN CALIBRATION: SER vs SNR for impulse channel, no FFE

n_symbols = 300_000
PAM = 4

# Generate PAM4 symbols: {-3, -1, +1, +3}
symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
# Normalize to {-0.75, -0.25, +0.25, +0.75}
symbols = symbols_raw * 0.25

sample_symb_fx = arrayFixedInt(X_W, X_F, symbols)
sample_symb_fx_arr = np.array([s.fValue for s in sample_symb_fx])
sym = np.round(sample_symb_fx_arr /0.25).astype(int)

# Impulse channel
x = symbols.copy()

snr_dbs = np.arange(0, 22, 2)
ser_sim = []
ber_sim = []

rng = np.random.default_rng(0)

print("Simulación de SER/BER para canal impulso + ruido AWGN, sin FFE. " \
"Quiero comparar con la curva teórica")
bits_true = np.array([gray_bits[v] for v in sym], dtype=np.uint8)

for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(x**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(x)) * np.sqrt(noise_var)

    y = x + noise
    #quantization 
    y = arrayFixedInt(X_W, X_F, y)

    ydec = np.array([slicer_fx(i, thr1=0.50, lvl1=0.25, lvl3=0.75) for i in y])
    ydec_slic = np.array([i.fValue for i in ydec])
    ydec_int = np.round(ydec_slic / 0.25).astype(int)

    bits_hat  = np.array([gray_bits[v] for v in ydec_int], dtype=np.uint8)

    ber = np.mean(bits_true != bits_hat)
    ser = np.mean(ydec_slic != sample_symb_fx_arr)
    print(f"at {snr_db}dB: SER = {ser:.2e}, BER = {ber:.2e}")
    ser_sim.append(ser)
    ber_sim.append(ber)

ser_sim = np.array(ser_sim)
ber_sim = np.array(ber_sim)    

snr_lin = 10**(snr_dbs/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)


# plt.figure()
# plt.semilogy(snr_dbs, ser_sim, 'o-', label="Sim (impulse + AWGN)")
# plt.semilogy(snr_dbs, ser_theory, '-', label="Theory PAM4 AWGN")
# plt.grid(True, which="both")
# plt.xlabel("SNR (dB)")
# plt.ylabel("SER")
# plt.title("PAM4 SER vs SNR (Impulse channel, no LMS)")
# plt.legend()
# plt.show()

#b=2*(M−1)/(M*log2​(M))*​Q(sqrt(6*log2(M)Eb/(M^2-1)/N0)
#b=3/4*Qfun(sqrt())
g = np.sqrt(0.4 * snr_lin)
ber_theory = (3*Qfun(g) + 2*Qfun(3*g) - Qfun(5*g)) / 4

import matplotlib.ticker as mticker

fig, ax = plt.subplots()
ax.semilogy(snr_dbs, ber_sim, 'o-', label="Sim (impulse + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim, 'o-', label="Sim (impulse + AWGN): SER")
ax.semilogy(snr_dbs, ber_theory, '--', label="Theory PAM4 AWGN: BER")
ax.semilogy(snr_dbs, ser_theory, '-', label="Theory PAM4 AWGN: SER")

# --- minor ticks on log Y ---
ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=12))
ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=tuple(range(2,10)), numticks=100))
ax.yaxis.set_minor_formatter(mticker.NullFormatter())  # hide minor tick labels (cleaner)

# grid
ax.grid(True, which='major', color='k', linestyle='-')
ax.grid(True, which='minor', color='k', linestyle='-', alpha=0.2)

ax.set_xlabel("SNR (dB)")
ax.set_ylabel("Error Rate")
ax.set_title("PAM4 Error Rate vs SNR (Impulse channel, no LMS)")
ax.legend()
plt.show()

# %% AWGN CALIBRATION: SER/BER vs SNR for impulse channel, FFE ON

n_symbols = 300_000
PAM = 4

# Generate PAM4 symbols: {-3, -1, +1, +3}
symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
# Normalize to {-0.75, -0.25, +0.25, +0.75}
symbols = symbols_raw * 0.25

sample_symb_fx = arrayFixedInt(X_W, X_F, symbols)
sample_symb_fx_arr = np.array([s.fValue for s in sample_symb_fx])
sym = np.round(sample_symb_fx_arr /0.25).astype(int)

bits_true = np.array([gray_bits[v] for v in sym], dtype=np.uint8)

# Impulse channel
x = symbols.copy()

snr_dbs = np.arange(0, 22, 2)
ser_sim_fp = []
ber_sim_fp = []

rng = np.random.default_rng(0)

print("Simulación de SER/BER para canal impulso + ruido AWGN, con FFE. " \
"Quiero comparar con la curva teórica")
STARTUP_DELAY = 3*len(FFE)
CMA_COUNT = int(1e4)
for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(x**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(x)) * np.sqrt(noise_var)

    y = x + noise
    #quantization 
    y = arrayFixedInt(X_W, X_F, y)

    # FFE Parameters
    FFE_LEN = 21
    CENTRAL_TAP = FFE_LEN//2
    FFE = np.zeros(FFE_LEN)
    FFE[CENTRAL_TAP] = 1.0
    mu_ffe = 1e-3
    mu_cma = 1e-4

    # quantization of variables
    mu_fx = Q(mu_ffe, MU_W, MU_F)
    mu_cma = Q(mu_cma, MU_W, MU_F)

    FFE = arrayFixedInt(W_W, W_F, FFE)
    mem_in_data = np.zeros(FFE_LEN).tolist()
    mem_in_data = arrayFixedInt(X_W, X_F, mem_in_data)
    ydec = []

    for i, sample in enumerate(y):
        # Shift memory
        for k in range(FFE_LEN-1, 0, -1):
            mem_in_data[k].assign(mem_in_data[k-1])
        mem_in_data[0].assign(sample)

        # FFE output
        #out_ffe = sample
        out_ffe         = FIR_fx(mem_in_data,FFE)
        out_ffe         = Q(out_ffe, X_W, X_F)
        
        # decider
        #out_slicer = out_ffe
        out_slicer = slicer_fx(out_ffe, thr1=0.50, lvl1=0.25, lvl3=0.75)
        ydec.append(out_slicer)
        

        error_slicer = Q(out_ffe-out_slicer, ACC_W, ACC_F)
        if (i%50000==0):
            print(f"iteration {i}/{n_symbols}")
        if i > STARTUP_DELAY:
            if i < CMA_COUNT:
                FFE = CMA(FFE,mem_in_data, out_ffe, mu_cma)
            else:
                FFE = LMS_fx(FFE,mem_in_data, error_slicer,mu_fx)

    ydec_slic = np.array([y.fValue for y in ydec])
    ydec_int = np.round(ydec_slic /0.25).astype(int)
    bits_hat  = np.array([gray_bits[v] for v in ydec_int], dtype=np.uint8)

    start = int(50e3)
    ser = np.mean(ydec_slic[start+CENTRAL_TAP:] != sample_symb_fx_arr[start:n_symbols-CENTRAL_TAP])
    ber = np.mean(bits_hat[start+CENTRAL_TAP:] != bits_true[start:n_symbols-CENTRAL_TAP])
    print(f"at {snr_db} dB: SER={ser:.2e}, BER={ber:.2e}")
    ser_sim_fp.append(ser)
    ber_sim_fp.append(ber)

ser_sim_fp = np.array(ser_sim_fp)
ber_sim_fp = np.array(ber_sim_fp)

snr_dbs_th = np.arange(0,22,0.1)
snr_lin = 10**(np.array(snr_dbs_th)/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)

g = np.sqrt(0.4 * snr_lin)
ber_theory = (3*Qfun(g) + 2*Qfun(3*g) - Qfun(5*g)) / 4

import matplotlib.ticker as mticker

fig, ax = plt.subplots()
ax.semilogy(snr_dbs, ber_sim_fp, 'o-', label="Sim (FFE + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim_fp, 'o-', label="Sim (FFE + AWGN): SER")
ax.semilogy(snr_dbs_th, ber_theory, '--', label="Theory PAM4 AWGN: BER")
ax.semilogy(snr_dbs_th, ser_theory, '-', label="Theory PAM4 AWGN: SER")

# --- minor ticks on log Y ---
ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=12))
ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=tuple(range(2,10)), numticks=100))
ax.yaxis.set_minor_formatter(mticker.NullFormatter())  # hide minor tick labels (cleaner)

# grid
ax.grid(True, which='major', color='k', linestyle='-')
ax.grid(True, which='minor', color='k', linestyle='-', alpha=0.2)

ax.set_xlabel("SNR (dB)")
ax.set_ylabel("Error Rate")
ax.set_title("PAM4 Error Rate vs SNR (Impulse channel, LMS)")
ax.legend()
plt.show()

# %% FFE VERIFICATION: SER/BER vs SNR for Low Pass channel, AWGN

n_symbols = 600_000
PAM = 4

Es_pam = (PAM**2 - 1) / 3  # 5

# Generate PAM4 symbols: {-3, -1, +1, +3}
symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
# Normalize to {-0.75, -0.25, +0.25, +0.75}
symbols = symbols_raw * 0.25


b,delay_ch = channel_fir(fcut=BW, fs_ch=SR, plt_en=False)
channel_symbols = np.convolve(b, symbols, mode="full")
channel_symbols = channel_symbols[delay_ch: delay_ch + len(symbols)]
n_samples = len(channel_symbols)

sample_symb_fx = arrayFixedInt(X_W, X_F, channel_symbols)
sample_symb_fx_arr = np.array([s.fValue for s in sample_symb_fx])

symb_fx = arrayFixedInt(X_W, X_F, symbols)
symb_fx_arr = np.array([s.fValue for s in symb_fx])

x = channel_symbols.copy()

snr_dbs = np.arange(6, 22, 2)
#snr_dbs = [20]
ser_sim_fp = []
ser_sim_ch = []
ber_sim_fp = []
ber_sim_ch = []

rng = np.random.default_rng(0)

print("Simulación de SER/BER para canal con att de 36db en 2GHz + ruido AWGN, con FFE. " \
"Quiero comparar con la curva teórica")
STARTUP_DELAY = 3*len(FFE)
CMA_COUNT = int(200e3)
for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(x**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(x)) * np.sqrt(noise_var)

    y = x + noise
    #quantization 
    y = arrayFixedInt(X_W, X_F, y)

    # FFE Parameters
    FFE_LEN = 21
    CENTRAL_TAP = FFE_LEN//2
    FFE = np.zeros(FFE_LEN)
    FFE[CENTRAL_TAP] = 1.0
    mu_ffe = 1e-3
    mu_cma = 1e-4

    # quantization of variables
    mu_fx = Q(mu_ffe, MU_W, MU_F)
    mu_cma = Q(mu_cma, MU_W, MU_F)

    FFE = arrayFixedInt(W_W, W_F, FFE)
    mem_in_data = np.zeros(FFE_LEN).tolist()
    mem_in_data = arrayFixedInt(X_W, X_F, mem_in_data)
    ydec = []
    ydec_ch = []
    FFE_history = []
    out_ffe_scope = []

    for i, sample in enumerate(y):
        if (i%50000==0):
            print(f"iteration {i}/{n_samples}")
        # Shift memory
        for k in range(FFE_LEN-1, 0, -1):
            mem_in_data[k].assign(mem_in_data[k-1])
        mem_in_data[0].assign(sample)

        # FFE output
        #out_ffe = sample
        out_ffe         = FIR_fx(mem_in_data,FFE)
        out_ffe         = Q(out_ffe, X_W, X_F)
        out_ffe_scope.append(out_ffe)
        # decider
        #out_slicer = out_ffe
        out_slicer = slicer_fx(out_ffe, thr1=0.50, lvl1=0.25, lvl3=0.75)
        ydec.append(out_slicer)

        #slicer channel
        slicer_ch = slicer_fx(sample, thr1=0.50, lvl1=0.25, lvl3=0.75)
        ydec_ch.append(slicer_ch)

        error_slicer = Q(out_ffe-out_slicer, ACC_W, ACC_F)
        if i > STARTUP_DELAY:
            if i < CMA_COUNT:
                FFE = CMA(FFE,mem_in_data, out_ffe, mu_cma)
            else:
                FFE = LMS_fx(FFE,mem_in_data, error_slicer,mu_fx)
            FFE_history.append(np.array([s.fValue for s in FFE]))

    ydec_slic = np.array([y.fValue for y in ydec])
    ydec_slic_ch = np.array([y.fValue for y in ydec_ch])
    out_ffe_scope = np.array([s.fValue for s in out_ffe_scope])
    
    ydec_int = np.round(ydec_slic /0.25).astype(int)
    ydec_int_ch = np.round(ydec_slic_ch /0.25).astype(int)
    bits_hat  = np.array([gray_bits[v] for v in ydec_int], dtype=np.uint8)
    bits_hat_ch  = np.array([gray_bits[v] for v in ydec_int_ch], dtype=np.uint8)

    start = int(75e3)
    ser = np.mean(ydec_slic[start+CENTRAL_TAP:] != symb_fx_arr[start:n_symbols-CENTRAL_TAP])
    ber = np.mean(bits_hat[start+CENTRAL_TAP:] != bits_true[start:n_symbols-CENTRAL_TAP])
    print(f"at {snr_db} dB: SER={ser:.2e}, BER={ber:.2e}")
    ser_sim_fp.append(ser)
    ber_sim_fp.append(ber)

    ser_ch = np.mean(ydec_slic_ch[start:] != symb_fx_arr[start:n_symbols])
    ber_ch = np.mean(bits_hat_ch[start:] != bits_true[start:n_symbols])
    print(f"at {snr_db} dB, Channel + AWGN only: SER={ser_ch:.2e}, BER={ber_ch:.2e}")
    ser_sim_ch.append(ser_ch)
    ber_sim_ch.append(ber_ch)

ser_sim_fp = np.array(ser_sim_fp)
ser_sim_ch = np.array(ser_sim_ch)

ber_sim_fp = np.array(ber_sim_fp)
ber_sim_ch= np.array(ber_sim_ch)


snr_dbs_th = np.arange(5,20,0.1)
snr_lin = 10**(np.array(snr_dbs_th)/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)

g = np.sqrt(0.4 * snr_lin)
ber_theory = (3*Qfun(g) + 2*Qfun(3*g) - Qfun(5*g)) / 4

import matplotlib.ticker as mticker

# SER/BER vs SNR plot
fig, ax = plt.subplots()
ax.semilogy(snr_dbs, ber_sim_fp, 'o-', label="Sim (FFE + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim_fp, 'o-', label="Sim (FFE + AWGN): SER")
ax.semilogy(snr_dbs, ber_sim_ch, '*-', label="Sim (Channel + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim_ch, '*-', label="Sim (Channel + AWGN): SER")
ax.semilogy(snr_dbs_th, ber_theory, '--', label="Theory PAM4 AWGN: BER")
ax.semilogy(snr_dbs_th, ser_theory, '-', label="Theory PAM4 AWGN: SER")

# --- minor ticks on log Y ---
ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=12))
ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=tuple(range(2,10)), numticks=100))
ax.yaxis.set_minor_formatter(mticker.NullFormatter())  # hide minor tick labels (cleaner)

# grid
ax.grid(True, which='major', color='k', linestyle='-')
ax.grid(True, which='minor', color='k', linestyle='-', alpha=0.2)

ax.set_xlabel("SNR (dB)")
ax.set_ylabel("Error Rate")
ax.set_title("PAM4 Error Rate vs SNR (LP channel, LMS)")
ax.legend()
plt.show()

# Bode plot of Heq, FFE, Channel
FFE_float = np.array([s.fValue for s in FFE])
heq = np.convolve(FFE_float, b, "full")

w, H_ch = sig.freqz(b, 1, worN=4096, fs=SR)
w, H_ffe = sig.freqz(FFE_float, 1, worN=4096, fs=SR)
w, H_eq = sig.freqz(heq, 1, worN=4096, fs=SR)


plt.figure(figsize=(10, 6))
plt.plot(w/1e9, 20*np.log10(np.abs(H_ch) + 1e-300))
plt.plot(w/1e9, 20*np.log10(np.abs(H_ffe) + 1e-300))
plt.plot(w/1e9, 20*np.log10(np.abs(H_eq) + 1e-300))
plt.title("Bode Plot of Butterworth Channel")
plt.ylabel("Magnitude [dB]")
plt.xlabel("Frequency [GHz]")
plt.grid(True)
plt.show()

# FFE vs time
plot_ffe(FFE_history, 10, title="FFE Coeff. Evolution - LP Channel")

# %% FFE VERIFICATION: SER/BER vs SNR for Butterworth channel, AWGN
FFE_LEN = 21
CENTRAL_TAP = FFE_LEN//2

n_symbols = 600_000
PAM = 4

Es_pam = (PAM**2 - 1) / 3  # 5
# Generate PAM4 symbols: {-3, -1, +1, +3}
symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
# Normalize to {-0.75, -0.25, +0.25, +0.75}
symbols = symbols_raw * 0.25


b,a,delay_ch = channel_butter(1.995e9, 1.85e9, SR, plt_en=True)
channel_symbols = sig.lfilter(b, a, symbols)
channel_symbols = channel_symbols[delay_ch: delay_ch + len(symbols)]
n_samples = len(channel_symbols)

sample_symb_fx = arrayFixedInt(X_W, X_F, channel_symbols)
sample_symb_fx_arr = np.array([s.fValue for s in sample_symb_fx])

symb_fx = arrayFixedInt(X_W, X_F, symbols)
symb_fx_arr = np.array([s.fValue for s in symb_fx])

x = channel_symbols.copy()

snr_dbs = np.arange(6, 22, 2)
#snr_dbs = [20]
ser_sim_fp = []
ser_sim_ch = []
ber_sim_fp = []
ber_sim_ch = []

rng = np.random.default_rng(0)

print("Simulación de SER/BER para canal con att de 53db en 2GHz + ruido AWGN, con FFE. " \
"Quiero comparar con la curva teórica")
CMA_COUNT = int(200e3)
for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(x**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(x)) * np.sqrt(noise_var)

    y = x + noise
    #quantization 
    y = arrayFixedInt(X_W, X_F, y)

    # FFE Parameters
    FFE_LEN = 21
    CENTRAL_TAP = FFE_LEN//2
    FFE = np.zeros(FFE_LEN)
    FFE[CENTRAL_TAP] = 1.0

    mu_ffe = 1e-3
    mu_cma = 1e-3

    STARTUP_DELAY = 5*len(FFE)
    # quantization of variables
    mu_fx = Q(mu_ffe, MU_W, MU_F)
    mu_cma = Q(mu_cma, MU_W, MU_F)

    FFE = arrayFixedInt(W_W, W_F, FFE)
    mem_in_data = np.zeros(FFE_LEN).tolist()
    mem_in_data = arrayFixedInt(X_W, X_F, mem_in_data)
    ydec = []
    ydec_ch = []
    FFE_history = []
    out_ffe_scope = []

    for i, sample in enumerate(y):
        if (i%50000==0):
            print(f"iteration {i}/{n_samples}")
        # Shift memory
        for k in range(FFE_LEN-1, 0, -1):
            mem_in_data[k].assign(mem_in_data[k-1])
        mem_in_data[0].assign(sample)

        # FFE output
        #out_ffe = sample
        out_ffe         = FIR_fx(mem_in_data,FFE)
        out_ffe         = Q(out_ffe, X_W, X_F)
        out_ffe_scope.append(out_ffe)
        # decider
        #out_slicer = out_ffe
        out_slicer = slicer_fx(out_ffe, thr1=0.50, lvl1=0.25, lvl3=0.75)
        ydec.append(out_slicer)

        #slicer channel
        slicer_ch = slicer_fx(sample, thr1=0.50, lvl1=0.25, lvl3=0.75)
        ydec_ch.append(slicer_ch)

        error_slicer = Q(out_ffe-out_slicer, ACC_W, ACC_F)
        if i > STARTUP_DELAY:
            if i < CMA_COUNT:
                FFE = CMA(FFE,mem_in_data, out_ffe, mu_cma)
            else:
                FFE = LMS_fx(FFE,mem_in_data, error_slicer,mu_fx)
            FFE_history.append(np.array([s.fValue for s in FFE]))

    ydec_slic = np.array([y.fValue for y in ydec])
    ydec_slic_ch = np.array([y.fValue for y in ydec_ch])
    out_ffe_scope = np.array([s.fValue for s in out_ffe_scope])
    
    ydec_int = np.round(ydec_slic /0.25).astype(int)
    ydec_int_ch = np.round(ydec_slic_ch /0.25).astype(int)
    bits_hat  = np.array([gray_bits[v] for v in ydec_int], dtype=np.uint8)
    bits_hat_ch  = np.array([gray_bits[v] for v in ydec_int_ch], dtype=np.uint8)

    start = int(75e3)
    ser = np.mean(ydec_slic[start+CENTRAL_TAP:] != symb_fx_arr[start:n_symbols-CENTRAL_TAP])
    ber = np.mean(bits_hat[start+CENTRAL_TAP:] != bits_true[start:n_symbols-CENTRAL_TAP])
    print(f"at {snr_db} dB: SER={ser:.2e}, BER={ber:.2e}")
    ser_sim_fp.append(ser)
    ber_sim_fp.append(ber)

    ser_ch = np.mean(ydec_slic_ch[start:] != symb_fx_arr[start:n_symbols])
    ber_ch = np.mean(bits_hat_ch[start:] != bits_true[start:n_symbols])
    print(f"at {snr_db} dB, Channel + AWGN only: SER={ser_ch:.2e}, BER={ber_ch:.2e}")
    ser_sim_ch.append(ser_ch)
    ber_sim_ch.append(ber_ch)

ser_sim_fp = np.array(ser_sim_fp)
ser_sim_ch = np.array(ser_sim_ch)

ber_sim_fp = np.array(ber_sim_fp)
ber_sim_ch= np.array(ber_sim_ch)


snr_dbs_th = np.arange(0,22,0.1)
snr_lin = 10**(np.array(snr_dbs_th)/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)

g = np.sqrt(0.4 * snr_lin)
ber_theory = (3*Qfun(g) + 2*Qfun(3*g) - Qfun(5*g)) / 4

import matplotlib.ticker as mticker

# SER/BER vs SNR Plot
fig, ax = plt.subplots()
ax.semilogy(snr_dbs, ber_sim_fp, 'o-', label="Sim (FFE + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim_fp, 'o-', label="Sim (FFE + AWGN): SER")
ax.semilogy(snr_dbs, ber_sim_ch, '*-', label="Sim (Channel + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim_ch, '*-', label="Sim (Channel + AWGN): SER")
ax.semilogy(snr_dbs_th, ber_theory, '--', label="Theory PAM4 AWGN: BER")
ax.semilogy(snr_dbs_th, ser_theory, '-', label="Theory PAM4 AWGN: SER")

# --- minor ticks on log Y ---
ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=12))
ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=tuple(range(2,10)), numticks=100))
ax.yaxis.set_minor_formatter(mticker.NullFormatter())  # hide minor tick labels (cleaner)

# grid
ax.grid(True, which='major', color='k', linestyle='-')
ax.grid(True, which='minor', color='k', linestyle='-', alpha=0.2)

ax.set_xlabel("SNR (dB)")
ax.set_ylabel("Error Rate")
ax.set_title("PAM4 Error Rate vs SNR (Butterworth channel, LMS)")
ax.legend()
plt.show()

# Bode plot of Heq, FFE, Channel
FFE_float = np.array([s.fValue for s in FFE])
heq = sig.lfilter(b, a, FFE_float)

w, H_ch = sig.freqz(b, a, worN=4096, fs=SR)
w, H_ffe = sig.freqz(FFE_float, 1, worN=4096, fs=SR)
w, H_eq = sig.freqz(heq, 1, worN=4096, fs=SR)

plt.figure(figsize=(10, 6))
plt.plot(w/1e9, 20*np.log10(np.abs(H_ch) + 1e-300))
plt.plot(w/1e9, 20*np.log10(np.abs(H_ffe) + 1e-300))
plt.plot(w/1e9, 20*np.log10(np.abs(H_eq) + 1e-300))
plt.title("Bode Plot of Butterworth Channel")
plt.ylabel("Magnitude [dB]")
plt.xlabel("Frequency [GHz]")
plt.grid(True)
plt.show()

# FFE vs time
plot_ffe(FFE_history, 10, title="FFE Coeff. Evolution - Butterworth Channel")
plt.show()
# %% FFE VERIFICATION: SER/BER vs SNR for another FIR channel, AWGN
FFE_LEN = 21
CENTRAL_TAP = FFE_LEN//2

n_symbols = 600_000
PAM = 4


# # Generate PAM4 symbols: {-3, -1, +1, +3}
# symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
# # Normalize to {-0.75, -0.25, +0.25, +0.75}
# symbols = symbols_raw * 0.25

# # APPLY CHANNEl
# b = channel_fir_nyquist_loss( nyq_loss_db=114, NTAPS=11, plt_en=True)
# a = 1

# channel_symbols = np.convolve(symbols, b, mode="full")
# channel_symbols = channel_symbols[:len(symbols)]
# n_samples = len(channel_symbols)

# sample_symb_fx = arrayFixedInt(X_W, X_F, channel_symbols)
# sample_symb_fx_arr = np.array([s.fValue for s in sample_symb_fx])

# symb_fx = arrayFixedInt(X_W, X_F, symbols)
# symb_fx_arr = np.array([s.fValue for s in symb_fx])

x = channel_symbols.copy()

#snr_dbs = np.arange(6, 22, 2)
snr_dbs = [20]
ser_sim_fp = []
ser_sim_ch = []
ber_sim_fp = []
ber_sim_ch = []

rng = np.random.default_rng(0)

print("Simulación de SER/BER para canal con att de 12db en 2GHz + ruido AWGN, con FFE. " \
"Quiero comparar con la curva teórica")
CMA_COUNT = int(200e3)
# SIMULATION
for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(x**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(x)) * np.sqrt(noise_var)

    y = x + noise
    #quantization 
    y = arrayFixedInt(X_W, X_F, y)

    # FFE Parameters
    FFE_LEN = 21
    CENTRAL_TAP = FFE_LEN//2
    FFE = np.zeros(FFE_LEN)
    FFE[CENTRAL_TAP] = 1.0

    mu_ffe = 1e-3
    mu_cma = 1e-3

    STARTUP_DELAY = 5*len(FFE)
    # quantization of variables
    mu_fx = Q(mu_ffe, MU_W, MU_F)
    mu_cma = Q(mu_cma, MU_W, MU_F)

    FFE = arrayFixedInt(W_W, W_F, FFE)
    mem_in_data = np.zeros(FFE_LEN).tolist()
    mem_in_data = arrayFixedInt(X_W, X_F, mem_in_data)
    ydec = []
    ydec_ch = []
    FFE_history = []
    out_ffe_scope = []
    error_scope = []

    for i, sample in enumerate(y):
        if (i%50000==0):
            print(f"iteration {i}/{n_samples}")
        # Shift memory
        for k in range(FFE_LEN-1, 0, -1):
            mem_in_data[k].assign(mem_in_data[k-1])
        mem_in_data[0].assign(sample)

        # FFE output
        #out_ffe = sample
        out_ffe         = FIR_fx(mem_in_data,FFE, ACC_W, ACC_F)
        out_ffe         = Q(out_ffe, X_W, X_F)
        out_ffe_scope.append(out_ffe)
        # decider
        #out_slicer = out_ffe
        out_slicer = slicer_fx(out_ffe, thr1=0.50, lvl1=0.25, lvl3=0.75)
        ydec.append(out_slicer)

        #slicer channel
        slicer_ch = slicer_fx(sample, thr1=0.50, lvl1=0.25, lvl3=0.75)
        ydec_ch.append(slicer_ch)

        error_slicer = Q(out_ffe-out_slicer, ACC_W, ACC_F)
        error_scope.append(error_slicer)    
        if i > STARTUP_DELAY:
            if i < CMA_COUNT:
                FFE = CMA(FFE,mem_in_data, out_ffe, mu_cma, W_W, W_F)
            else:
                FFE = LMS_fx(FFE,mem_in_data, error_slicer,mu_fx, W_W, W_F)
            FFE_history.append(np.array([s.fValue for s in FFE]))

    ydec_slic = np.array([y.fValue for y in ydec])
    ydec_slic_ch = np.array([y.fValue for y in ydec_ch])
    out_ffe_scope = np.array([s.fValue for s in out_ffe_scope])
    error_scope = np.array([s.fValue for s in error_scope])

    ydec_int = np.round(ydec_slic /0.25).astype(int)
    ydec_int_ch = np.round(ydec_slic_ch /0.25).astype(int)
    bits_hat  = np.array([gray_bits[v] for v in ydec_int], dtype=np.uint8)
    bits_hat_ch  = np.array([gray_bits[v] for v in ydec_int_ch], dtype=np.uint8)

    start = int(200e3)
    ser = np.mean(ydec_slic[start+CENTRAL_TAP:] != symb_fx_arr[start:n_symbols-CENTRAL_TAP])
    ber = np.mean(bits_hat[start+CENTRAL_TAP:] != bits_true[start:n_symbols-CENTRAL_TAP])
    print(f"at {snr_db} dB: SER={ser:.2e}, BER={ber:.2e}")
    ser_sim_fp.append(ser)
    ber_sim_fp.append(ber)

    ser_ch = np.mean(ydec_slic_ch[start:] != symb_fx_arr[start:n_symbols])
    ber_ch = np.mean(bits_hat_ch[start:] != bits_true[start:n_symbols])
    print(f"at {snr_db} dB, Channel + AWGN only: SER={ser_ch:.2e}, BER={ber_ch:.2e}")
    ser_sim_ch.append(ser_ch)
    ber_sim_ch.append(ber_ch)

ser_sim_fp = np.array(ser_sim_fp)
ser_sim_ch = np.array(ser_sim_ch)

ber_sim_fp = np.array(ber_sim_fp)
ber_sim_ch= np.array(ber_sim_ch)


snr_dbs_th = np.arange(0,22,0.1)
snr_lin = 10**(np.array(snr_dbs_th)/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)

g = np.sqrt(0.4 * snr_lin)
ber_theory = (3*Qfun(g) + 2*Qfun(3*g) - Qfun(5*g)) / 4

import matplotlib.ticker as mticker

# SER/BER vs SNR Plot
fig, ax = plt.subplots()
ax.semilogy(snr_dbs, ber_sim_fp, 'o-', label="Sim (FFE + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim_fp, 'o-', label="Sim (FFE + AWGN): SER")
ax.semilogy(snr_dbs, ber_sim_ch, '*-', label="Sim (Channel + AWGN): BER")
ax.semilogy(snr_dbs, ser_sim_ch, '*-', label="Sim (Channel + AWGN): SER")
ax.semilogy(snr_dbs_th, ber_theory, '--', label="Theory PAM4 AWGN: BER")
ax.semilogy(snr_dbs_th, ser_theory, '-', label="Theory PAM4 AWGN: SER")

# --- minor ticks on log Y ---
ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=12))
ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=tuple(range(2,10)), numticks=100))
ax.yaxis.set_minor_formatter(mticker.NullFormatter())  # hide minor tick labels (cleaner)

# grid
ax.grid(True, which='major', color='k', linestyle='-')
ax.grid(True, which='minor', color='k', linestyle='-', alpha=0.2)

ax.set_xlabel("SNR (dB)")
ax.set_ylabel("Error Rate")
ax.set_title("PAM4 Error Rate vs SNR (FIR channel, LMS)")
ax.legend()
plt.show()

# Bode plot of Heq, FFE, Channel
FFE_float = np.array([s.fValue for s in FFE])
heq = sig.lfilter(b, a, FFE_float)

w, H_ch = sig.freqz(b, a, worN=4096, fs=SR)
w, H_ffe = sig.freqz(FFE_float, 1, worN=4096, fs=SR)
w, H_eq = sig.freqz(heq, 1, worN=4096, fs=SR)

plt.figure(figsize=(10, 6))
plt.plot(w/1e9, 20*np.log10(np.abs(H_ch) + 1e-300))
plt.plot(w/1e9, 20*np.log10(np.abs(H_ffe) + 1e-300))
plt.plot(w/1e9, 20*np.log10(np.abs(H_eq) + 1e-300))
plt.title("Bode Plot of FIR Channel")
plt.ylabel("Magnitude [dB]")
plt.xlabel("Frequency [GHz]")
plt.grid(True)
plt.show()

# FFE vs time
plot_ffe(FFE_history, 10, title="FFE Coeff. Evolution - FIR Channel")
plt.show()

# %%
# Plot 1: Input vs Output scatter
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# Transmitted symbols
downsample=10
Nplot = min(5000, n_symbols)
axs[0, 0].plot(symb_fx_arr[::downsample], linestyle='None', marker='.', markersize=4, alpha=0.7)
axs[0, 0].axhline(0.75, color='r', linestyle='--', alpha=0.5)
axs[0, 0].axhline(0.25, color='r', linestyle='--', alpha=0.5)
axs[0, 0].axhline(-0.25, color='r', linestyle='--', alpha=0.5)
axs[0, 0].axhline(-0.75, color='r', linestyle='--', alpha=0.5)
axs[0, 0].set_title("Transmitted Symbols (PAM4)")
axs[0, 0].set_ylabel("Amplitude")
axs[0, 0].grid(True, alpha=0.3)

# Received signal (quantized)
y_float = [s.fValue for s in y]
axs[0, 1].plot(y_float[::downsample], linestyle='None', marker='.', markersize=4, alpha=0.7)
axs[0, 1].axhline(0.5, color='g', linestyle=':', alpha=0.5)
axs[0, 1].axhline(0.0, color='g', linestyle=':', alpha=0.5)
axs[0, 1].axhline(-0.5, color='g', linestyle=':', alpha=0.5)
axs[0, 1].set_title("Received Signal (Q8.7 Quantized)")
axs[0, 1].set_ylabel("Amplitude")
axs[0, 1].grid(True, alpha=0.3)

# FFE output
axs[1, 0].plot(out_ffe_scope[::downsample], linestyle='None', marker='.', markersize=4, alpha=0.7)
axs[1, 0].axhline(0.75, color='r', linestyle='--', alpha=0.5)
axs[1, 0].axhline(0.25, color='r', linestyle='--', alpha=0.5)
axs[1, 0].axhline(-0.25, color='r', linestyle='--', alpha=0.5)
axs[1, 0].axhline(-0.75, color='r', linestyle='--', alpha=0.5)
axs[1, 0].set_title("FFE Output (Equalized)")
axs[1, 0].set_ylabel("Amplitude")
axs[1, 0].set_xlabel("Sample Index")
axs[1, 0].grid(True, alpha=0.3)

# Error
axs[1, 1].plot(error_scope[::downsample], linewidth=0.5, alpha=0.7)
axs[1, 1].set_title("Slicer Error")
axs[1, 1].set_ylabel("Error")
axs[1, 1].set_xlabel("Sample Index")
axs[1, 1].grid(True, alpha=0.3)

## FFE HISTORY
fig1 = plot_ffe(FFE_history, 1000, f"FFE Coefficient Evolution - mu={mu_ffe:.1e}")

# get last FFE taps
FFE = FFE_history[-1]
fig, ax = plt.subplots(figsize=(10, 5))
ax.stem(range(FFE_LEN), FFE, basefmt=' ')
ax.set_title("Final FFE Tap Weights")
ax.set_xlabel("Tap Index")
ax.set_ylabel("Coefficient Value")
ax.axvline(CENTRAL_TAP, color='r', linestyle='--', alpha=0.5, label='Center Tap')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()

# Bode plot of Heq, FFE, Channel
heq = sig.lfilter(b, 1, FFE)

w, H_ch = sig.freqz(b, a, worN=4096, fs=SR)
w, H_ffe = sig.freqz(FFE, 1, worN=4096, fs=SR)
w, H_eq = sig.freqz(heq, 1, worN=4096, fs=SR)

plt.figure(figsize=(10, 6))
plt.plot(w/1e9, 20*np.log10(np.abs(H_ch) + 1e-300), label="Channel")
plt.plot(w/1e9, 20*np.log10(np.abs(H_ffe) + 1e-300), label="FFE")
plt.plot(w/1e9, 20*np.log10(np.abs(H_eq) + 1e-300), label="H equivalent")
plt.title("Bode Plot of Channel")
plt.ylabel("Magnitude [dB]")
plt.xlabel("Frequency [GHz]")
plt.legend()
plt.grid(True)
plt.show()

#%% 
n_symbols = 600_000
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


#%%
FFE_LEN = 21 
FFE = np.array([-6.19411469e-04, -1.46055222e-03, -2.39253044e-04,  8.14795494e-04,
       -4.23955917e-03,  5.16009331e-03, -1.77979469e-03,  1.05035305e-03,
       -4.50968742e-04,  7.29918480e-03,  2.48351741e+00, -1.21298206e+00,
       -2.70409226e-01, -1.16955400e-01, -6.51751757e-02, -4.35215235e-02,
        3.49032402e-01, -6.46786690e-02, -4.25608158e-02, -2.47280598e-02,
       -1.85259581e-02])

#%%
x = channel_symbols.copy()
x = arrayFixedInt(X_W, X_F, x)
# FFE Parameters
#%%
FFE_LEN = 21
CENTRAL_TAP = FFE_LEN//2
FFE = np.zeros(FFE_LEN)
FFE[CENTRAL_TAP] = 8388607/2**(W_F)

mu_ffe = 1e-3
mu_cma = 1/2**10

# quantization of variables
mu_fx = Q(mu_ffe, MU_W, MU_F)
mu_cma = Q(mu_cma, MU_W, MU_F)

FFE_fx = arrayFixedInt(W_W, W_F, FFE)


# simulation
mem_in_data = np.zeros(FFE_LEN).tolist()
mem_in_data = arrayFixedInt(X_W, X_F, mem_in_data)
mem_in_data_list = []
ydec = []
out_ffe_scope = []
cma_error_scope = []
FFE_history = []
ffe_updates = []
updates = []
R = DeFixedInt(W_W, W_F, roundMode='trunc', saturateMode='saturate')
aux = np.array([-0.75, -0.25, 0.25, 0.75])
# R.value = 1.64
R.value = np.mean(np.abs(aux)**4) / np.mean(np.abs(aux)**2)

CMA_COUNT = 700_000
STARTUP_DELAY = 3*len(FFE)-1
for i, sample in enumerate(x):
    if (i%50000==0):
        print(f"iteration {i}/{n_samples}")
    # Shift memory
    for k in range(FFE_LEN-1, 0, -1):
        mem_in_data[k].assign(mem_in_data[k-1])
    mem_in_data[0].assign(sample)

    mem_in_data_list.append([s.value for s in mem_in_data])

    if (i==STARTUP_DELAY):
        print(mem_in_data)
        mem_in_data_first = [s.fValue for s in mem_in_data]
        
    # FFE output
    #out_ffe = sample
    out_ffe         = FIR_fx(mem_in_data,FFE_fx, ACC_W, ACC_F)
    out_ffe         = Q(out_ffe, X_W, X_F)
    out_ffe_scope.append(out_ffe)
    # decider

    out_slicer = slicer_fx(out_ffe, thr1=0.50, lvl1=0.25, lvl3=0.75)
    ydec.append(out_slicer)

    # error calculation for CMA
    error_fx = Q(out_ffe*out_ffe - R, W_W, W_F)
    cma_error_scope.append(error_fx)

    # ffe_fx = np.zeros(FFE_LEN, dtype=object)
    # ffe_fx = arrayFixedInt(W_W, W_F, ffe_fx)
    # for k in range(len(FFE_fx)):
    #     upd = mem_in_data[k] * error_fx * mu_fx * out_ffe
    #     ffe_fx[k] = Q(FFE_fx[k] - upd, W_W, W_F)
    #     print(f"tap {k}: mem={mem_in_data[k]}, upd={upd}, new_tap={FFE_fx[k] - upd}")
    #     if (k == CENTRAL_TAP):
    #         updates.append(upd)
    # ffe_updates.append(ffe_fx)

    

    # print(f"iteration {i}")
    if i >= STARTUP_DELAY:
        if i < CMA_COUNT:
            FFE = CMA(FFE_fx,mem_in_data, out_ffe, mu_cma, W_W, W_F)
        else:
            FFE = LMS_fx(FFE_fx,mem_in_data, error_slicer,mu_fx)
    FFE_history.append(np.array([s for s in FFE_fx]))



#%%
bin_symb_ch = [str(s.value) for s in sample_symb_fx]
with open("./output_mem/channel_symbols_int.mem", "w") as f:
    for s in bin_symb_ch:
        f.write(s + "\n")

bin_symb_ch = [s.bit() for s in sample_symb_fx]
with open("./output_mem/channel_symbols.mem", "w") as f:
    for s in bin_symb_ch:
        f.write(s + "\n")


bin_symb = [str(s.value) for s in symb_fx]
with open("./output_mem/symbols_int.mem", "w") as f:
    for s in bin_symb:
        f.write(s + "\n")


bin_out_ffe = [str(s.value) for s in out_ffe_scope]
with open("./output_mem/out_ffe_int.mem", "w") as f:
    for s in bin_out_ffe:
        f.write(s + "\n")

bin_ydec = [str(s.value) for s in ydec]
with open("./output_mem/dec_symb_int.mem", "w") as f:
    for s in bin_ydec:
        f.write(s + "\n")

bin_cma_error = [str(s.value) for s in cma_error_scope]
with open("./output_mem/cma_error_int.mem", "w") as f:
    for s in bin_cma_error:
        f.write(s + "\n")

# %% READ FROM VERILOG SIMULATION OUTPUT
out_ffe_rtl = np.loadtxt("C:\\Users\\denis\\Documents\\beca\\porcom-tp-final\\verilog\\testbench\\fir_out_rtl.txt", dtype=int)
out_ffe_rtl = out_ffe_rtl[2:]

out_ffe_golden = np.array([s.value for s in out_ffe_scope])
out_ffe_golden = out_ffe_golden[0:len(out_ffe_rtl)]

print("Max abs diff:", np.max(np.abs(out_ffe_rtl - out_ffe_golden)))
# %%
# %% READ FROM VERILOG SIMULATION OUTPUT
out_slicer_rtl = np.loadtxt("C:\\Users\\denis\\Documents\\beca\\porcom-tp-final\\verilog\\testbench\\slicer_out_rtl.txt", dtype=int)
out_slicer_rtl = out_slicer_rtl[2:]

out_slicer_golden = np.array([s.value for s in ydec])
out_slicer_golden = out_slicer_golden[0:len(out_slicer_rtl)]

print("Max abs diff:", np.max(np.abs(out_slicer_rtl - out_slicer_golden)))

# %%
# %% READ FROM VERILOG SIMULATION OUTPUT
cma_error_rtl = np.loadtxt("C:\\Users\\denis\\Documents\\beca\\porcom-tp-final\\verilog\\testbench\\cma_error_rtl.txt", dtype=int)
cma_error_rtl = cma_error_rtl[2:]

cma_error_golden = np.array([s.value for s in cma_error_scope])
cma_error_golden = cma_error_golden[0:len(cma_error_rtl)]

print("Max abs diff:", np.max(np.abs(cma_error_rtl - cma_error_golden)))

# %%
# %% READ FROM VERILOG SIMULATION OUTPUT
out_ffe_rtl = np.loadtxt("C:\\Users\\denis\\Documents\\beca\\porcom-tp-final\\verilog\\testbench\\out_ffe.txt", dtype=int)
out_ffe_rtl = out_ffe_rtl[1:]

out_ffe_golden = np.array([s.value for s in out_ffe_scope])
if out_ffe_rtl.shape[0] < out_ffe_golden.shape[0]:
    out_ffe_golden = out_ffe_golden[0:len(out_ffe_rtl)]
else:
    out_ffe_rtl = out_ffe_rtl[0:len(out_ffe_golden)]

print("Max abs diff:", np.max(np.abs(out_ffe_rtl - out_ffe_golden)))
for i in range(len(out_ffe_rtl)):
    if out_ffe_rtl[i] != out_ffe_golden[i]:
        print(f"Mismatch at index {i}: RTL={out_ffe_rtl[i]}, Golden={out_ffe_golden[i]}")
# %%
