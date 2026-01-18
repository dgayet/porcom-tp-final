import plotly.graph_objs as go
from plotly.subplots import make_subplots
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import control as cm
import scipy.fft as fft
import scipy.signal as sig
from scipy.signal import welch
from scipy import signal
from scipy.signal import lfilter as filt
from numpy.random import randint as randi
import scipy.signal as scF
from tool.rcosdesign import rcosdesign as rcos
from matplotlib.gridspec import GridSpec
import matplotlib
from fractions import Fraction
from decimal import *
from tool._fixedInt import *
from scipy.stats import kurtosis, skew
from numpy import i0
matplotlib.use('TkAgg')

def slicer(x: float, M: int = 4):
    """
    Slices a received voltage 'x' to the nearest PAM-M level.
    Assumes M is an even integer >= 2.
    The ideal levels are: -(M-1), ..., -3, -1, +1, +3, ..., (M-1).

    Args:
        x: The received analog voltage level (float).
        M: The number of modulation levels (e.g., 4 for PAM-4).
    Returns:
        The sliced integer symbol level.
    """
    # 1. Validate M
    if M < 2 or M % 2 != 0:
        raise ValueError("M must be an even integer greater than or equal to 2.")

    # 2. Find the nearest odd integer
    # We find the nearest integer, then check if it's odd or even.
    # Python 3's round() rounds .5 to the nearest *even* number, 
    # e.g., round(0.5)=0, round(1.5)=2. This logic correctly 
    # handles the decision boundaries (which are at 0, +/-2, +/-4, ...).
    
    r = int(np.round(x))
    
    if r % 2 != 0:
        # If r is already odd (e.g., 1, 3, -1), it's our value.
        sliced_val = r
    else:
        # If r is even (e.g., 0, 2, -2), it's a decision boundary.
        if x > r:
            sliced_val = r + 1
        else:
            sliced_val = r - 1

    # 3. Clamp the result to the valid PAM-M range
    max_level = M - 1
    # Clamp to the minimum level: -(M-1)
    output = max(sliced_val, -max_level)
    # Clamp to the maximum level: (M-1)
    output = min(output, max_level)
    
    return output

def GET_SER(slicer_scope, CENTRAL_TAP, symbols, start_ber=200000):
    # SYMBOL ERROR RATE
    slicer_arr = np.array(slicer_scope)
    align = CENTRAL_TAP
    if align != 0:
        slicer_align = slicer_arr[align:]
        symbols_trunc = symbols[:-align]
    else:
        slicer_align = slicer_arr
        symbols_trunc = symbols
    print(slicer_align.shape, symbols_trunc.shape)
    #start_ber = 200000
    if len(slicer_align) <= start_ber:
        # simulation not long enough
        print("Warning: not enough samples for BER start index. Lower start_ber or increase n_symbols.")
    else:
        errors_bool = slicer_align[start_ber:] != symbols_trunc[start_ber:]
        n_errors = np.sum(errors_bool)
        ncomp = len(slicer_align[start_ber:])
        SER_value = n_errors / ncomp
        print(f"SER: {SER_value:.6e}  ({n_errors}/{ncomp})")
    
    return SER_value

def FIR(samples,coeffs):
    return np.dot(samples,coeffs)

def LMS(fir,samples,error,mu):
    fir = fir - samples*error*mu
    return fir

def CMA(fir,samples,yk,mu, R=1.64):
    error=yk**2 - R
    fir = fir - samples*error*yk*mu
    return fir

def MCMA_PAM4(fir, samples, yk, mu, R_levels=(0.2, 1.8)):
    y2 = yk*yk
    # choose closest target energy (inner vs outer)
    R = R_levels[0] if abs(y2 - R_levels[0]) < abs(y2 - R_levels[1]) else R_levels[1]
    error = y2 - R
    fir = fir - samples * error * yk * mu
    return fir

def db(data):
    return 20*np.log10(np.abs(data))

def db10(data):
    return 10*np.log10(np.abs(data))

def channel_fir(fcut, fs_ch, plt_en=False):

    print("CHANNEL DESIGN: FIR")

    ORDER = 128  # filter length - 1
    f_cut = 0.48 * fs_ch  # Hz
    b = sig.firwin(ORDER+1, f_cut/(fs_ch/2), window='blackmanharris')
    a = 1  # FIR → no denominator

    # Frequency response
    w, h = sig.freqz(b, a, worN=4096, fs=fs_ch)
    f_target = fcut

    # Attenuation at Nyquist frequency
    idx = np.argmin(np.abs(w - f_target))
    attenGHz = 20 * np.log10(np.abs(h[idx]))
    print(f"Attenuation at {f_target/1e9:.2f} GHz: {attenGHz:.2f} dB\n")


    # PLOT CHANNEL FREQUENCY RESPONSE
    if plt_en:
        mag_db = 20*np.log10(np.abs(h) + 1e-300)          # avoid log(0)
        phase = np.unwrap(np.angle(h))

        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax[0].plot(w/1e9, mag_db)
        ax[0].set_title("Bode Plot of FIR Channel Filter")
        ax[0].set_ylabel("Magnitude [dB]")
        ax[0].set_ylim(-100, 20)
        ax[0].grid(True)

        ax[1].plot(w/1e9, phase)
        ax[1].set_xlabel("Frequency [GHz]")
        ax[1].set_ylabel("Phase [rad]")
        ax[1].grid(True)

        ax[0].set_xlim(0, 2.1)   # up to 10 GHz (adjust to your SR)
        plt.show()

    delay = ORDER // 2

    return b, delay

def channel_butter(fsb, fp, fs_ch, plt_en=False):
    """
    Butterworth lowpass channel.
    fcut: cutoff in Hz
    fs_ch: sampling frequency in Hz
    order: IIR order (typical 3..8)
    Returns: (b, a, delay_est)
    """
    print("CHANNEL DESIGN: Butterworth IIR")

    #Fs   = 4e9
    #fp   = 1.85e9     # passband edge (choose based on what you want to preserve)
    gpass = 0.5      # dB ripple allowed in passband (Butterworth is monotonic, this is "max loss")
    gstop = 20       # dB attenuation required in stopband

    N, Wn = sig.buttord(wp=fp, ws=fsb, gpass=gpass, gstop=gstop, fs=fs_ch)
    b, a  = sig.butter(N, Wn, btype="low", fs=fs_ch)

    print("Butter order =", N, "Wn =", Wn)

    w, h = sig.freqz(b, a, worN=16384, fs=fs_ch)

    for f in [1.9e9, 1.95e9, 1.99e9, 2.0e9]:
        idx = np.argmin(np.abs(w - f))
        att = 20*np.log10(np.abs(h[idx]) + 1e-300)
        print(f"Attenuation at {w[idx]/1e9:.3f} GHz: {att:.2f} dB")

    if plt_en:
        mag_db = 20*np.log10(np.abs(h) + 1e-300)          # avoid log(0)
        phase = np.unwrap(np.angle(h))

        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax[0].plot(w/1e9, mag_db)
        ax[0].set_title("Bode Plot of IIR Channel Filter")
        ax[0].set_ylabel("Magnitude [dB]")
        ax[0].set_ylim(-100, 20)
        ax[0].grid(True)

        ax[1].plot(w/1e9, phase)
        ax[1].set_xlabel("Frequency [GHz]")
        ax[1].set_ylabel("Phase [rad]")
        ax[1].grid(True)

        ax[0].set_xlim(0, 2.1)   # up to 10 GHz (adjust to your SR)
        plt.show()

    # Butterworth has frequency-dependent group delay.
    # A simple "best" sample alignment proxy is the peak of the impulse response.
    imp = np.zeros(1024)
    imp[0] = 1.0
    h_imp = sig.lfilter(b, a, imp)
    delay_est = int(np.argmax(np.abs(h_imp)))

    return b, a, delay_est


def channel_cheby(fcut, fs_ch, plt_en=False):
    """

    Chebyshev bandpass channel.
    fcut: cutoff in Hz
    fs_ch: sampling frequency in Hz
    order: IIR order (typical 3..8)
    Returns: (b, a, delay_est)
    """
    print("CHANNEL DESIGN: Butterworth IIR")

    # Digital cutoff (normalized to Nyquist)
    fp  = 1.7e9 / (fs_ch / 2)           # passband edge (choose based on what you want to preserve)
    fsb = (0.99*fcut)/ (fs_ch / 2)      # stopband edge (keep < 2 GHz; push closer if you want)
    rp  = 0.5        # dB passband ripple
    rs  = 36         # dB stopband attenuation

    N, Wn = sig.cheb1ord(wp=fp, ws=fsb, gpass=rp, gstop=rs, fs=fs_ch)
    b, a  = sig.cheby1(N, Wn, rp, btype='low', fs=fs_ch)
    # Frequency response
    w, h = sig.freqz(b, a, worN=4096, fs=fs_ch)

    # Attenuation at fcut (or any f_target you want)
    idx = np.argmin(np.abs(w - fcut))
    atten = 20 * np.log10(np.abs(h[idx]) + 1e-300)
    print(f"Attenuation at {fcut/1e9:.2f} GHz: {atten:.2f} dB\n")

    if plt_en:
        mag_db = 20*np.log10(np.abs(h) + 1e-300)          # avoid log(0)
        phase = np.unwrap(np.angle(h))

        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax[0].plot(w/1e9, mag_db)
        ax[0].set_title("Bode Plot of FIR Channel Filter")
        ax[0].set_ylabel("Magnitude [dB]")
        ax[0].set_ylim(-100, 20)
        ax[0].grid(True)

        ax[1].plot(w/1e9, phase)
        ax[1].set_xlabel("Frequency [GHz]")
        ax[1].set_ylabel("Phase [rad]")
        ax[1].grid(True)

        ax[0].set_xlim(0, 2.1)   # up to 10 GHz (adjust to your SR)
        plt.show()

    # Butterworth has frequency-dependent group delay.
    # A simple "best" sample alignment proxy is the peak of the impulse response.
    imp = np.zeros(1024)
    imp[0] = 1.0
    h_imp = sig.lfilter(b, a, imp)
    delay_est = int(np.argmax(np.abs(h_imp)))

    return b, a, delay_est


def rrc(CHANNEL_UP, symbols, n_symbols):
    rcos_filt = rcos(0.1, 40, CHANNEL_UP, 1,'sqrt')
    # uncomment if DSPtools is being used
    # (_,rcos_filt) = rcos(0.1, 1, CHANNEL_UP, 40,'sqrt')

    rcos_delay = len(rcos_filt)//2 # Filter delay

    symbols_up = np.zeros(n_symbols*CHANNEL_UP) # Oversampled symbol sequence
    symbols_up[::CHANNEL_UP] = symbols

    out_upsampled = np.convolve(symbols_up, rcos_filt, mode="full")
    return out_upsampled[rcos_delay : rcos_delay + len(symbols_up)]   # valid samples

def impulse_response(b, FFE, CHANNEL_UP, plt_en):

    # Create large delta
    delta_k = np.zeros(int(10e3))
    delta_k[0] = 1.0

    # Channel impulse response at high sampling rate
    model_fir = sig.lfilter(b, 1, delta_k)

    # Upsample FFE to channel rate
    resampled_FFE = scF.resample_poly(FFE, CHANNEL_UP, 1)
    resampled_FFE /= np.sum(resampled_FFE)

    # Compute combined channel + FFE response
    combined_response_FFE = sig.lfilter(b, 1, resampled_FFE)

    # Find timing phase (integer offset)
    phase = int(np.argmax(np.abs(combined_response_FFE)) % CHANNEL_UP)

    if plt_en:
        fig4, (ax00,ax01,ax02) = plt.subplots(3, 1, figsize=(10,8))

        n = np.arange(len(FFE))

        # FFE taps
        ax00.stem(n, FFE, basefmt=" ", linefmt='b', markerfmt='bo', label="FFE taps")
        ax00.grid(True)
        ax00.legend()

        # Combined channel + FFE response (downsampled)
        comb = combined_response_FFE[phase::CHANNEL_UP][:len(FFE)]
        ax01.stem(n, comb, basefmt=" ", linefmt='r', markerfmt='ro', label="Combined Response")
        ax01.grid(True)
        ax01.legend()

        # Channel-only response (downsampled at same phase!)
        model_down = model_fir[phase::CHANNEL_UP][:len(FFE)]
        ax02.stem(n, model_down, basefmt=" ", linefmt='g', markerfmt='go', label="Channel Only")
        ax02.grid(True)
        ax02.legend()

        plt.tight_layout()
        plt.show()

    return phase

def plot_error(error_scope, DOWN_PLOT,title=""):
    fig = plt.figure()
    plt.plot(error_scope[::DOWN_PLOT])
    plt.title(title)
    plt.grid(True)
    plt.show()
    return fig

def plot_symbols(in_ffe_scope, ffe_out_scope, DOWN_PLOT, title):
    fig5, (ax50,ax51) = plt.subplots(2, 1)

    markerline51, stemlines51, baseline51 = ax51.stem(ffe_out_scope[::DOWN_PLOT],'r',label="Out FFE")

    plt.setp(stemlines51, visible=False)   # no vertical lines
    baseline51.set_visible(False)          # no baseline
    ax51.grid(True)
    ax51.legend()

    markerline50, stemlines50, baseline50 = ax50.stem(in_ffe_scope[::DOWN_PLOT],'b',label="In FFE")
    plt.setp(stemlines50, visible=False)   # no vertical lines
    baseline50.set_visible(False)          # no baseline
    ax50.grid(True)
    ax50.legend()
    plt.title(title)
    
    fig5.tight_layout()
    plt.show()
    return fig5

def plot_ffe_frec(FFE, fs, fig=None):
    w, h = sig.freqz(FFE, worN=4096, fs=fs)  # w in Hz
    if fig==None:
        #plt.figure(figsize=(10,6))
        fig, ax = plt.subplots(figsize=(10,6))
    else:
        ax = fig.gca()
    ax.plot(w/1e9, 20*np.log10(np.abs(h)))
    ax.set_title("Bode Plot of FFE Rx")
    ax.set_ylabel("Magnitude [dB]")
    ax.set_xlabel("frequency [GHz]")
    # ax.set_xlim(0, 10)   # show only up to 100 GHz
    # ax.set_ylim(-100, 0)   # show only up to 100 GHz
    ax.grid(True)
    #ax.show()
    return fig

def plot_ffe(FFE_history, DOWN_PLOT, title="FFE Coefficient Evolution"):
    fig = plt.figure()
    plt.plot(FFE_history[::DOWN_PLOT])
    plt.title(title)
    plt.xlabel("Iteration")
    plt.ylabel("Tap value")
    plt.grid(True)
    #plt.show()
    return fig

