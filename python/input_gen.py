#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.special import erfc
# from fxpmath import Fxp
import math

################
#### DEFINITIONS
################

def slicer_pam4(x: float):
    """
    PAM4 slicer for symbols at {-0.75, -0.25, +0.25, +0.75}
    with decision thresholds at {-0.5, 0, +0.5}

    Args:
        x: The received analog voltage level (float).
    Returns:
        The sliced symbol level.
    """
    if x < -0.5:
        return -0.75
    elif x < 0:
        return -0.25
    elif x < 0.5:
        return 0.25
    else:
        return 0.75

def GET_SER(slicer_scope, CENTRAL_TAP, symbols, start_ber=200):
    """Calculate Symbol Error Rate"""
    slicer_arr = np.array(slicer_scope)
    align = CENTRAL_TAP
    slicer_align = slicer_arr[align:]
    symbols_trunc = symbols[:-align]

    if len(slicer_align) <= start_ber:
        print(f"Warning: not enough samples for SER calculation (need >{start_ber})")
        return 0.0

    errors_bool = slicer_align[start_ber:] != symbols_trunc[start_ber:]
    n_errors = np.sum(errors_bool)
    ncomp = len(slicer_align[start_ber:])
    SER_value = n_errors / ncomp
    print(f"SER: {SER_value:.6e}  ({n_errors}/{ncomp})")

    return SER_value

def FIR(samples, coeffs):
    """FIR filter: sum of products"""
    return np.dot(samples, coeffs)

def LMS(fir, samples, error, mu):
    """LMS adaptation: subtract gradient"""
    fir = fir - samples * error * mu
    return fir

def CMA(fir,samples,yk,mu):
    # R=8.2/(4**2)
    a = np.array([-0.75, -0.25, 0.25, 0.75])
    R = np.mean(np.abs(a)**4) / np.mean(np.abs(a)**2)
    error=yk**2-R
    fir = fir - samples*error*yk*mu
    return fir

def channel_fir(fcut, fs_ch, plt_en=False):
    """Create FIR channel model"""
    print("CHANNEL DESIGN: FIR")

    ORDER = 128  # filter length - 1
    f_cut = 0.49 * fs_ch  # Hz
    b = sig.firwin(ORDER+1, f_cut/(fs_ch/2), window='blackmanharris')
    a = 1  # FIR → no denominator

    # Frequency response
    w, h = sig.freqz(b, a, worN=4096, fs=fs_ch)
    f_target = fcut

    # Attenuation at Nyquist frequency
    idx = np.argmin(np.abs(w - f_target))
    attenGHz = 20 * np.log10(np.abs(h[idx]))
    print(f"Attenuation at {f_target/1e9:.2f} GHz: {attenGHz:.2f} dB\n")

    if plt_en:
        plt.figure(figsize=(10,6))
        plt.plot(w/1e9, 20*np.log10(np.abs(h)))
        plt.title("Bode Plot of FIR Channel Filter")
        plt.xlabel("Frequency [GHz]")
        plt.ylabel("Magnitude [dB]")
        plt.xlim(0, 10)
        plt.ylim(-100, 0)
        plt.grid(True)
        plt.show()

    delay = ORDER // 2
    return b, delay

def channel_fir_nyquist_loss_dc_loss(
    nyq_loss_db=20,
    NTAPS=31,
    plt_en=False
):

    w = np.linspace(0, np.pi, 8192)
    f_norm = w / np.pi

    mag_db = -nyq_loss_db * (f_norm ** 0.5)
    mag = 10**(mag_db / 20)

    H = mag * np.exp(1j * 0)

    h = np.fft.irfft(H)
    h = h[:NTAPS]

    h = sig.minimum_phase(h, method='homomorphic')

    # --- Force Nyquist attenuation ---
    w_plot, H_plot = sig.freqz(h, worN=4096)
    nyq_idx = np.argmin(np.abs(w_plot - np.pi))
    scale = 10**(-nyq_loss_db/20) / np.abs(H_plot[nyq_idx])
    h *= scale

    nyq_att_db = 20*np.log10(np.abs(H_plot[nyq_idx]) * scale)
    print(f"[CHANNEL] Attenuation at Nyquist: {nyq_att_db:.2f} dB")

    if plt_en:
        plt.plot(w_plot/np.pi, 20*np.log10(np.abs(H_plot)*scale))
        plt.axvline(1.0, color='r', linestyle='--')
        plt.axhline(-nyq_loss_db, color='g', linestyle='--')
        plt.grid()
        plt.show()

        plt.stem(h)
        plt.grid()
        plt.show()

    return h

def channel_fir_nyquist_loss(
    nyq_loss_db=20,
    NTAPS=31,
    plt_en=False
):
    """
    Baud-rate discrete-time channel model
    - 1 sample / symbol
    - minimum-phase
    - DC gain normalized to 1
    - ~nyq_loss_db attenuation at Nyquist
    """

    # ---------- Frequency grid (0 → Nyquist) ----------
    w = np.linspace(0, np.pi, 8192)
    f_norm = w / np.pi   # 0 → 1

    # ---------- Loss model (smooth, monotonic) ----------
    mag_db = -nyq_loss_db * (f_norm ** 0.5)
    mag = 10**(mag_db / 20)

    # Zero-phase spectrum
    H = mag.astype(np.complex128)

    # ---------- Impulse response ----------
    h = np.fft.irfft(H)
    h = h[:NTAPS]

    # Minimum-phase equivalent
    h = sig.minimum_phase(h, method='homomorphic')

    # ---------- Normalize DC gain ----------
    h /= np.sum(h)

    # ---------- Log DC and Nyquist ----------
    w_plot, H_plot = sig.freqz(h, worN=4096)

    dc_gain_db  = 20*np.log10(np.abs(H_plot[0]))
    nyq_idx     = np.argmin(np.abs(w_plot - np.pi))
    nyq_att_db  = 20*np.log10(np.abs(H_plot[nyq_idx]))

    print(f"[CHANNEL] DC gain: {dc_gain_db:.2f} dB")
    print(f"[CHANNEL] Attenuation @ Nyquist (f_norm=1): {nyq_att_db:.2f} dB")

    # ---------- Optional plots ----------
    if plt_en:
        plt.figure(figsize=(8,4))
        plt.plot(w_plot/np.pi, 20*np.log10(np.abs(H_plot)))
        plt.axvline(1.0, color='r', linestyle='--', label="Nyquist")
        plt.axhline(0, color='k', linestyle='--', label="DC")
        plt.xlabel("Normalized Frequency (f / Nyquist)")
        plt.ylabel("Magnitude [dB]")
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(8,4))
        try:
            plt.stem(h, use_line_collection=True)
        except TypeError:
            plt.stem(h)
        plt.title("Baud-rate channel impulse response")
        plt.xlabel("Symbol index")
        plt.ylabel("Amplitude")
        plt.grid()
        plt.tight_layout()
        plt.show()

    return h

def SER(tx, rx, discard=200):
    tx = np.asarray(tx)
    rx = np.asarray(rx)

    # discard startup / transients
    tx = tx[discard:]
    rx = rx[discard:]

    n = min(len(tx), len(rx))
    tx = tx[:n]
    rx = rx[:n]

    n_errors = np.sum(tx != rx)
    ser = n_errors / n

    print(f"SER = {ser:.3e}  ({n_errors}/{n})")

    return ser, n_errors, n


def slicer_pam4_vec(x):
    """
    Vectorized PAM4 slicer
    """
    x = np.asarray(x)

    y = np.empty_like(x)
    y[x < -0.5]              = -0.75
    y[(x >= -0.5) & (x < 0)] = -0.25
    y[(x >= 0) & (x < 0.5)]  = 0.25
    y[x >= 0.5]              = 0.75

    return y


def theoretical_ser_pam4(snr_db):
    """
    Calculate theoretical SER for PAM4 in AWGN channel.

    For PAM4 with Gray coding, the standard formula is:
    SER = 2*(M-1)/M * Q(sqrt(6*SNR / (M²-1)))

    For M=4:
    SER = (3/2) * Q(sqrt(6*SNR/15))
        = (3/2) * Q(sqrt(0.4*SNR))

    where Q(x) = (1/2)*erfc(x/sqrt(2)) and SNR = Es/N0

    Args:
        snr_db: Signal-to-Noise Ratio in dB (Es/N0 basis)

    Returns:
        Theoretical Symbol Error Rate
    """
    # Convert SNR from dB to linear
    snr_linear = 10**(snr_db / 10)

    # For PAM-M with Gray coding, the standard formula is:
    # SER = 2*(M-1)/M * Q(sqrt(6*SNR / (M²-1)))
    #
    # For M=4:
    # SER = 2*3/4 * Q(sqrt(6*SNR/15))
    #     = (3/2) * Q(sqrt(0.4*SNR))
    #
    # where SNR = Es/N0 (symbol energy to noise PSD ratio)

    argument = np.sqrt(0.4 * snr_linear)

    # Q function: Q(x) = 0.5 * erfc(x / sqrt(2))
    q_func = 0.5 * erfc(argument / np.sqrt(2))

    # For PAM4 with Gray coding:
    ser_theoretical = (3/2) * q_func

    return ser_theoretical


################
#### SIMULATION
################

# Parameters

only_noise = 0 #set to 1 to bypass channel

n_symbols = int(3e7)
PAM = 4
SR = 28e9          # symbol rate (baud)
BR = SR * np.log2(PAM)   # bit rate
BW = SR/2  # Nyquist BandWidth

print(f"=" * 60)
print(f"STARTING SIMULATION: PAM {PAM} EQUALIZER")
print(f"=" * 60)
print(f"\tBit Rate: {BR/1e9:.2f} Gbps")
print(f"\tSymbol Rate: {SR/1e9:.2f} GBd")
print(f"\tBand Width: {BW/1e9:.2f} GHz\n")

# Generate PAM4 symbols: {-3, -1, +1, +3}
symbols_raw = 2*np.random.randint(0, PAM, n_symbols) - PAM + 1

# Normalize to {-0.75, -0.25, +0.25, +0.75}
symbols = symbols_raw * 0.25

################
#### SAVE FOR RTL
################

# Quantize to Q8.7 for RTL simulation
# sym_q87 = Fxp(
#     symbols,
#     signed=True,
#     n_word=8,
#     n_frac=7,
#     rounding='trunc',
#     overflow='saturate'
# )

# sym_int = sym_q87.val.astype(np.int16)

# output_path = "/home/ignaciobalbo/serdes_digital_ffe/modules/ffe_serial/model/fixed_point/sym_q87.txt"
# np.savetxt(output_path, sym_int, fmt="%d")
# print(f"Saved {len(sym_int)} samples to: {output_path}\n")

# continue symbol gen

print(f"Symbol levels: {np.unique(symbols)}")
print(f"Expected: [-0.75 -0.25  0.25  0.75]\n")

# Channel model
# b, delay_ch = channel_fir(fcut=BW, fs_ch=SR, plt_en=False)
b = channel_fir_nyquist_loss( nyq_loss_db=114,
                              NTAPS=11,
                              plt_en=False)
# print(f"Channel delay: {delay_ch} samples\n")

# Apply channel (bypassed when only_noise=1)
if only_noise:
    print("*** CHANNEL BYPASSED - ONLY NOISE APPLIED ***\n")
    channel_symbols = symbols
else:
    # channel_symbols = np.convolve(b, symbols, mode="full")
    # channel_symbols = channel_symbols[delay_ch : delay_ch + len(symbols)]
    channel_symbols = np.convolve(symbols, b, mode="full")
    channel_symbols = channel_symbols[:len(symbols)]

# Add noise
snr_db = 15
print(f"SNR: {snr_db} dB\n")

snr_lin = 10**(snr_db / 10)
Ps = np.mean(channel_symbols**2)
# Standard definition: SNR = Es/N0, where N0 = 2*σ² for real signals
# Therefore: σ² = Es/(2*SNR)
noise_var = Ps / (2 * snr_lin)
noise = np.random.randn(len(channel_symbols)) * np.sqrt(noise_var)

rx_channel = channel_symbols + noise

print(f"Signal power (Es): {Ps:.6f}")
print(f"Noise variance (σ²): {noise_var:.6f}")
print(f"N0 (= 2*σ²): {2*noise_var:.6f}")
print(f"Actual SNR (Es/N0): {10*np.log10(Ps/(2*noise_var)):.2f} dB\n")

################
#### SAVE FOR RTL
################

# Quantize to Q8.7 for RTL simulation
# rx_q87 = Fxp(
#     rx_channel,
#     signed=True,
#     n_word=8,
#     n_frac=7,
#     rounding='trunc',
#     overflow='saturate'
# )

# rx_int = rx_q87.val.astype(np.int16)

# output_path = "/home/ignaciobalbo/serdes_digital_ffe/modules/ffe_serial/model/fixed_point/rx_q87_dec.txt"
# np.savetxt(output_path, rx_int, fmt="%d")
# print(f"Saved {len(rx_int)} samples to: {output_path}\n")

################
#### THEORETICAL SER CALCULATION
################

print("=" * 60)
print("THEORETICAL SER CALCULATION")
print("=" * 60)

ser_theory = theoretical_ser_pam4(snr_db)

print(f"Theoretical SER (approximation): {ser_theory:.6e}")

################
#### DIRECT SER (NO EQUALIZATION)
################

print("=" * 60)
print("DIRECT SLICING (NO EQUALIZATION)")
print("=" * 60)

rx_sliced_direct = slicer_pam4_vec(rx_channel)
ser_direct, n_errors_direct, n_total = SER(symbols, rx_sliced_direct, discard=1000)
print()

################
#### FFE EQUALIZATION
################

print("=" * 60)
print("FFE EQUALIZATION")
print("=" * 60)

FFE_LEN = 21
FFE = np.zeros(FFE_LEN)
CENTRAL_TAP = FFE_LEN // 2
FFE[CENTRAL_TAP] = 1.0  # Initialize with unity at center tap

mu_ffe = 1e-3
mu_cma = mu_ffe
STARTUP_DELAY = 5 * FFE_LEN

mem_in_data = np.zeros(FFE_LEN)

error_scope = []
slicer_scope = []
ffe_out_scope = []
FFE_history = []

counter_adap = 0
PARALLELISM = 64
CMA_COUNT = 1.5e7

print("Running FFE adaptation...")
for ii, sample_data in enumerate(rx_channel):
    # Shift delay line
    mem_in_data[1:] = mem_in_data[:-1]
    mem_in_data[0] = sample_data

    # FIR filtering
    out_ffe = FIR(mem_in_data, FFE)

    # Slicer (decision device)
    out_slicer = slicer_pam4(out_ffe)

    # Calculate error
    error_slicer = out_ffe - out_slicer

    # Store for analysis
    error_scope.append(error_slicer)
    ffe_out_scope.append(out_ffe)
    slicer_scope.append(out_slicer)

    # LMS adaptation (after startup delay)
    if ii > STARTUP_DELAY:
        if counter_adap > PARALLELISM:
            if ii < CMA_COUNT:
                FFE = CMA(FFE,mem_in_data,out_ffe,mu_cma)
            else:# if counter_adap >= FFE_LEN:
                FFE = LMS(FFE, mem_in_data, error_slicer, mu_ffe)

            counter_adap = 0
        else:
            counter_adap+=1
        FFE_history.append(FFE.copy())

print("Adaptation complete!\n")

# Calculate SER with FFE
ser_ffe = GET_SER(slicer_scope, CENTRAL_TAP, symbols, start_ber=int(2e7))
print()

################
#### COMPARISON TABLE
################

print("=" * 60)
print("SER COMPARISON SUMMARY")
print("=" * 60)
print(f"{'Metric':<40} {'Value':>15}")
print("-" * 60)
print(f"{'SNR (dB)':<40} {snr_db:>15.2f}")
print(f"{'Number of symbols':<40} {n_symbols:>15,}")
print("-" * 60)
print(f"{'Theoretical SER (approximation)':<40} {ser_theory:>15.6e}")
print("-" * 60)
print(f"{'Simulated SER (no equalization)':<40} {ser_direct:>15.6e}")
print(f"{'  Errors':<40} {n_errors_direct:>15,}")
print(f"{'  Samples':<40} {n_total:>15,}")
print("-" * 60)
print(f"{'Simulated SER (with FFE)':<40} {ser_ffe:>15.6e}")
print("-" * 60)
print(f"{'Ratio (simulated/theoretical)':<40} {ser_direct/ser_theory:>15.2f}")
print("=" * 60)
print()

################
#### PLOTTING
################

# Print final FFE coefficients
print(f"Final FFE coefficients:")
print(FFE)
print()

# Plot 1: Input vs Output scatter
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# Transmitted symbols
downsample=100
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
axs[0, 1].plot(rx_q87[::downsample], linestyle='None', marker='.', markersize=4, alpha=0.7)
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

plt.tight_layout()
plt.savefig('/home/ignaciobalbo/serdes_digital_ffe/modules/ffe_serial/model/fixed_point/results_serial/pam4_signals.png', dpi=150, bbox_inches='tight')
print("Saved plot: pam4_signals.png")
plt.show()

# Plot 2: FFE coefficient evolution
if len(FFE_history) > 0:
    FFE_history_array = np.array(FFE_history)

    fig, ax = plt.subplots(figsize=(12, 6))
    for tap in range(FFE_LEN):
        # Downsample for visibility
        downsample = 1000
        ax.plot(FFE_history_array[::downsample, tap],
                label=f'Tap {tap}' if tap == CENTRAL_TAP else '',
                linewidth=2 if tap == CENTRAL_TAP else 0.5,
                alpha=1.0 if tap == CENTRAL_TAP else 0.5)

    ax.set_title("FFE Coefficient Adaptation")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Coefficient Value")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.savefig('/home/ignaciobalbo/serdes_digital_ffe/modules/ffe_serial/model/fixed_point/results_serial/ffe_adaptation.png', dpi=150, bbox_inches='tight')
    print("Saved plot: ffe_adaptation.png")
    plt.show()

# Plot 3: Final FFE taps
fig, ax = plt.subplots(figsize=(10, 5))
ax.stem(range(FFE_LEN), FFE, basefmt=' ')
ax.set_title("Final FFE Tap Weights")
ax.set_xlabel("Tap Index")
ax.set_ylabel("Coefficient Value")
ax.axvline(CENTRAL_TAP, color='r', linestyle='--', alpha=0.5, label='Center Tap')
ax.legend()
ax.grid(True, alpha=0.3)
plt.savefig('/home/ignaciobalbo/serdes_digital_ffe/modules/ffe_serial/model/fixed_point/results_serial/ffe_taps.png', dpi=150, bbox_inches='tight')
print("Saved plot: ffe_taps.png")
plt.show()

# Plot 4: SER vs SNR comparison
print("\nGenerating SER vs SNR curve...")
snr_range = np.linspace(5, 25, 21)
ser_theory_curve = [theoretical_ser_pam4(snr) for snr in snr_range]

fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogy(snr_range, ser_theory_curve, 'b-', linewidth=2, label='Theoretical (exact)')
ax.semilogy(snr_db, ser_theory, 'bo', markersize=10, label=f'Theoretical @ {snr_db} dB')
ax.semilogy(snr_db, ser_direct, 'rs', markersize=10, label=f'Simulated (no EQ) @ {snr_db} dB')
ax.semilogy(snr_db, ser_ffe, 'g^', markersize=10, label=f'Simulated (FFE) @ {snr_db} dB')

ax.set_xlabel('SNR (dB)')
ax.set_ylabel('Symbol Error Rate (SER)')
ax.set_title('PAM4 SER Performance: Theory vs Simulation')
ax.grid(True, which='both', alpha=0.3)
ax.legend()
ax.set_xlim(5, 25)
ax.set_ylim(1e-8, 1)
plt.savefig('/home/ignaciobalbo/serdes_digital_ffe/modules/ffe_serial/model/fixed_point/results_serial/ser_comparison.png', dpi=150, bbox_inches='tight')
print("Saved plot: ser_comparison.png")
plt.show()

print("\n" + "=" * 60)
print("SIMULATION COMPLETE")
print("=" * 60)
# %%
