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

# SYMBOL GENERATOR
n_symbols = int(1e6)

# MODULATION
PAM = 4
SR = 4e9          # symbol rate (baud)
BR = SR * np.log2(PAM)   # bit rate

BW = SR/2 # Nyquist BandWidth

print(f"- - - - - - - - - - - - - - - - - - - - -")
print(f"STARTING SIMULATION: PAM {PAM} EQUALIZER")
print(f"\tBit Rate: {BR/1e9:.2f}Gbps\n\tSymbol Rate: {SR/1e9:.2f}GBd\n\tBand Width: {BW/1e9:.2f}GHz\n")

# PAM-4 average symbol energy = (M^2 - 1) / 3 = 5
Es_pam = (PAM**2 - 1) / 3
norm = np.sqrt(Es_pam)   # normalizado a Es = 1

#%% BASE CASE: No channel, no noise, no pulse shaping -> FFE must be identity
# Symbol generation
symbols = (2*np.random.randint(0,PAM,n_symbols)-PAM+1)/norm

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

# Registers for FFE
mem = np.zeros(FFE_LEN)
ydec = np.zeros(n_symbols)

error_scope      = []
ffe_scope        = []
ffe_out_scope    = []
in_ffe_scope     = []
FFE_history     = []

for i, sample in enumerate(channel_symbols):
    # Shift memory
    mem[1:] = mem[:-1]
    mem[0] = sample

    # FFE output
    yffe = FIR(mem, FFE)

    # decider
    yslicer = slicer(yffe*norm, PAM)/norm
    ydec[i] = yslicer

    error_slicer = yffe - yslicer

    FFE = LMS(FFE, mem, error_slicer, mu_ffe)
    FFE_history.append(FFE.copy())

#snr=GET_SER(ydec, CENTRAL_TAP, symbols)
#display(snr)

print("Caso Base: el canal es la identidad, sin ruido, sin pulse shaping:" \
"espero que el FFE quede igual al impulso, solamente habrá que alinear la secuencia de entrada" \
"con la de salida del FFE (delay porque el tap central está en el medio del filtro)")

valid = np.arange(CENTRAL_TAP, n_symbols)
SER0 = np.mean(ydec[valid] != symbols[valid - CENTRAL_TAP])
print("SER base (impulso, sin ruido) =", SER0)
print("Central tap final =", FFE[CENTRAL_TAP])

fig1 = plot_ffe(FFE_history, 10, f"FFE Coefficient Evolution - mu={mu_ffe:.1e}")
# %% AWGN CALIBRATION: SER vs SNR for impulse channel, no FFE

# n_symbols = 1e7
n_symbols = 300_000
PAM = 4

Es_pam = (PAM**2 - 1) / 3  # 5
norm = np.sqrt(Es_pam)

# Es=1 symbols
symbols = (2*np.random.randint(0, PAM, n_symbols) - PAM + 1) / norm
symbols = symbols.astype(float)

# Impulse channel
x = symbols.copy()

snr_dbs = np.arange(-10, 20, 2)
ser_sim = []

rng = np.random.default_rng(0)

print("Simulación de SER para canal impulso + ruido AWGN, sin FFE. " \
"Quiero comparar con la curva teórica")

for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(x**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(x)) * np.sqrt(noise_var)

    y = x + noise

    ydec = np.array([slicer(v*norm, PAM)/norm for v in y])

    ser = np.mean(ydec != symbols)
    ser_sim.append(ser)

ser_sim = np.array(ser_sim)

snr_lin = 10**(snr_dbs/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)
plt.figure()
plt.semilogy(snr_dbs, ser_sim, 'o-', label="Sim (impulse + AWGN)")
plt.semilogy(snr_dbs, ser_theory, '-', label="Theory PAM4 AWGN")
plt.grid(True, which="both")
plt.xlabel("SNR (dB)")
plt.ylabel("SER")
plt.title("PAM4 SER vs SNR (Impulse channel, no LMS)")
plt.legend()
plt.show()

fig2 = plt.gcf()
#fig2.savefig("ser_pam4_impulse_awgn.png")
# %% AWGN CALIBRATION: SER vs SNR for impulse channel, FFE on

# FFE Parameters
FFE_LEN = 21
CENTRAL_TAP = FFE_LEN//2
mu_ffe = 1e-6

# n_symbols = 1e7
n_symbols = 300_000
PAM = 4

# Es=1 symbols
symbols = (2*np.random.randint(0, PAM, n_symbols) - PAM + 1) / norm
symbols = symbols.astype(float)

# Impulse channel
x = symbols.copy()

Es_pam = (PAM**2 - 1) / 3  # 5
norm = np.sqrt(Es_pam)


snr_dbs = np.arange(-10, 20, 2)
ser_sim = []

STARTUP_DELAY = 5*len(FFE)
CMA_COUNT = int(1e4)
for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(x**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(x)) * np.sqrt(noise_var)

    y = x + noise

    mem = np.zeros(FFE_LEN)
    ydec = np.zeros(n_symbols)

    FFE = np.zeros(FFE_LEN)
    FFE[CENTRAL_TAP] = 1.0
    FFE_history = []
    for i, sample in enumerate(y):
        # Shift memory
        mem[1:] = mem[:-1]
        mem[0] = sample

        # FFE output
        yffe = FIR(mem, FFE)

        # decider
        yslicer = slicer(yffe*norm, PAM)/norm
        ydec[i] = yslicer

        error_slicer = yffe - yslicer

        if (i > STARTUP_DELAY):
            if i < CMA_COUNT:
                FFE = CMA(FFE,mem,yffe,1e-4,1.64)
            else:
                FFE = LMS(FFE, mem, error_slicer, mu_ffe)
            FFE_history.append(FFE.copy())
    start = int(50e3)
    ser = np.mean(ydec[start+CENTRAL_TAP:] != symbols[start:n_symbols-CENTRAL_TAP])
    ser_sim.append(ser)

ser_sim = np.array(ser_sim)

snr_lin = 10**(snr_dbs/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)

plt.figure()
plt.semilogy(snr_dbs, ser_sim, 'o-', label="Sim (impulse + AWGN)")
plt.semilogy(snr_dbs, ser_theory, '-', label="Theory PAM4 AWGN")
plt.grid(True, which="both")
plt.xlabel("SNR (dB)")
plt.ylabel("SER")
plt.title("PAM4 SER vs SNR (Impulse channel, LMS)")
plt.legend()
plt.show()
# %% FFE Verification: Channel  + FFE
# No hay normalización, sino symbolos: 0.75, 0.25, -0.25, -0.75

n_symbols = int(1e6)
PAM = 4

# Generate PAM4 symbols: {-3, -1, +1, +3}
symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1
# Normalize to {-0.75, -0.25, +0.25, +0.75}
symbols = symbols_raw * 0.25


# APPLY CHANNEl
b = channel_fir_nyquist_loss( nyq_loss_db=114,
                              NTAPS=11,
                              plt_en=True)
# print(f"Channel delay: {delay_ch} samples\n")

channel_symbols = np.convolve(symbols, b, mode="full")
channel_symbols = channel_symbols[:len(symbols)]
#samples_symbols = channel_symbols/np.sqrt(np.mean(channel_symbols**2))# now sample at best integer point
#n_samples = len(samples_symbols)
n_samples = len(channel_symbols)
samples_symbols = channel_symbols.copy()

# FFE Parameters
FFE_LEN = 21
CENTRAL_TAP = FFE_LEN//2

mu_ffe = 1e-3
mu_cma = 1e-3
#%% NO AWGN
FFE = np.zeros(FFE_LEN)
FFE[CENTRAL_TAP] = 1.0

# Registers for FFE
mem = np.zeros(FFE_LEN)
ydec = np.zeros(n_samples)
ydec_ch_only = np.zeros(n_samples)

error_scope      = []
ffe_scope        = []
ffe_out_scope    = []
in_ffe_scope     = []
FFE_history     = []

STARTUP_DELAY = 5*len(FFE)
CMA_COUNT = 250000
for i, sample in enumerate(samples_symbols):
    # Shift memory
    mem[1:] = mem[:-1]
    mem[0] = sample

    # FFE output
    yffe = FIR(mem, FFE)
    ffe_out_scope.append(yffe)

    # decider
    yslicer = slicer_pam4_87(yffe)
    ydec[i] = yslicer
    ydec_ch_only[i] = slicer_pam4_87(sample)

    error_slicer = yffe - yslicer
    error_scope.append(error_slicer)
    if (i > STARTUP_DELAY):
        if i < CMA_COUNT:
            FFE = CMA(FFE,mem,yffe,mu_cma)
        else:
            FFE = LMS(FFE, mem, error_slicer, mu_ffe)
        FFE_history.append(FFE.copy())

    if i in [CMA_COUNT-1, CMA_COUNT, CMA_COUNT+1, CMA_COUNT+100]:
        print("i", i,
          "yffe_rms", np.sqrt(np.mean(np.array(ffe_out_scope[-2000:])**2)),
          "err_rms",  np.sqrt(np.mean(np.array(error_scope[-2000:])**2)),
          "FFE_center", FFE[CENTRAL_TAP])

print(GET_SER(ydec, CENTRAL_TAP, symbols))
print(GET_SER(ydec_ch_only, 0, symbols))

fig1 = plot_ffe(FFE_history, 1000, f"FFE Coefficient Evolution - mu={mu_ffe:.1e}")

# %% FFE Verification: Channel + AWGN + FFE
snr_dbs = np.array([20, 18, 15, 12, 10, 6])
ser_sim = []
ser_sim_ch = []
FFE_history = []

mu_ffe = 5e-2
mu_cma = 1e-3
STARTUP_DELAY = 5*len(FFE)
CMA_COUNT = 400_000
channel_symbols = samples_symbols.copy()

print("Simulación de SER para canal con ISI + ruido AWGN, con FFE LMS+CMA.")
print("--------------------- PARAMETROS ---------------------")
print(f"n_symbols: {n_symbols}")
print(f"mu_ffe: {mu_ffe}, mu_cma: {mu_cma}")
print(f"CMA_COUNT: {CMA_COUNT}, STARTUP_DELAY: {STARTUP_DELAY}")
print("------------------------------------------------------")


for snr_db in snr_dbs:
    snr_lin = 10**(snr_db/10)

    Ps = np.mean(channel_symbols**2)              
    noise_var = Ps / (2 * snr_lin)      
    noise = rng.standard_normal(len(channel_symbols)) * np.sqrt(noise_var)

    y = channel_symbols + noise

    mem = np.zeros(FFE_LEN)
    ydec = np.zeros(len(y))
    ydec_ch_only=np.zeros(len(y))

    FFE = np.zeros(FFE_LEN)
    FFE[CENTRAL_TAP] = 1.0
    for i, sample in enumerate(y):
        # Shift memory
        mem[1:] = mem[:-1]
        mem[0] = sample

        # FFE output
        yffe = FIR(mem, FFE)

        # decider
        yslicer = slicer_pam4_87(yffe)
        ydec[i] = yslicer
        ydec_ch_only[i] = slicer_pam4_87(y[i])

        error_slicer = yffe - yslicer
        if (i > STARTUP_DELAY):
            if i < CMA_COUNT:
                FFE = CMA(FFE,mem,yffe,mu_cma)
            else:
                FFE = LMS(FFE, mem, error_slicer, mu_ffe)

            if (snr_db == snr_dbs[0]):
                FFE_history.append(FFE.copy())

    print(f"ydec SER at SNR {snr_db} dB:")
    ser = GET_SER(ydec, CENTRAL_TAP, symbols, start_ber=400_000)
    ser_sim.append(ser)
    print(f"ydec_ch_only SER at SNR {snr_db} dB:")
    ser= GET_SER(ydec_ch_only, 0, symbols, start_ber=400_000)
    ser_sim_ch.append(ser)

ser_sim = np.array(ser_sim)
ser_sim_ch = np.array(ser_sim_ch)

snr_lin = 10**(snr_dbs/10)

ser_theory = 2 * (1 - 1/PAM) * Qfun(
    np.sqrt(6 * snr_lin / (PAM**2 - 1))
)
plt.figure()
plt.semilogy(snr_dbs, ser_sim, 'o-', label="Sim (Channel + FFE + AWGN)")
plt.semilogy(snr_dbs, ser_theory, '-', label="Theory PAM4 AWGN")
plt.semilogy(snr_dbs, ser_sim_ch, 'x--', label="Sim (Channel + AWGN)")
plt.grid(True, which="both")
plt.xlabel("SNR (dB)")
plt.ylabel("SER")
plt.title("PAM4 SER vs SNR (channel, LMS)")
plt.legend()
plt.show()

#%%
# Plot 1: Input vs Output scatter
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# Transmitted symbols
downsample=10
Nplot = min(5000, n_symbols)
axs[0, 0].plot(symbols[::downsample], linestyle='None', marker='.', markersize=4, alpha=0.7)
axs[0, 0].axhline(0.75, color='r', linestyle='--', alpha=0.5)
axs[0, 0].axhline(0.25, color='r', linestyle='--', alpha=0.5)
axs[0, 0].axhline(-0.25, color='r', linestyle='--', alpha=0.5)
axs[0, 0].axhline(-0.75, color='r', linestyle='--', alpha=0.5)
axs[0, 0].set_title("Transmitted Symbols (PAM4)")
axs[0, 0].set_ylabel("Amplitude")
axs[0, 0].grid(True, alpha=0.3)

# Received signal (quantized)
axs[0, 1].plot(samples_symbols[::downsample], linestyle='None', marker='.', markersize=4, alpha=0.7)
axs[0, 1].axhline(0.5, color='g', linestyle=':', alpha=0.5)
axs[0, 1].axhline(0.0, color='g', linestyle=':', alpha=0.5)
axs[0, 1].axhline(-0.5, color='g', linestyle=':', alpha=0.5)
axs[0, 1].set_title("Received Signal (Q8.7 Quantized)")
axs[0, 1].set_ylabel("Amplitude")
axs[0, 1].grid(True, alpha=0.3)

# FFE output
axs[1, 0].plot(ffe_out_scope[::downsample], linestyle='None', marker='.', markersize=4, alpha=0.7)
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

w, H_ch = sig.freqz(b, 1, worN=4096, fs=SR)
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


# %%
