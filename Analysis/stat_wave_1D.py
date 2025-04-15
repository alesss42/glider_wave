import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, detrend

def fast_psd(x, nfft_in, fs):
    """
    Compute a properly normalized power spectral density (PSD) for a time series.
    
    Parameters:
        x (np.array): Input data array.
        nfft_in (int): Desired FFT length.
        fs (float): Sampling frequency (Hz) for the PSD.
        
    Returns:
        powe (np.array): PSD values (one-sided, excluding the DC term).
        fre (np.array): Frequency vector corresponding to the PSD.
    """
    max_ind = len(x)
    nfft = min(nfft_in, max_ind)
    # Determine the number of segments (using at least 50% overlap as in the MATLAB code)
    repeats = int(2 * max_ind // nfft)
    if max_ind == nfft:
        repeats = 1
        if nfft != nfft_in:
            print(f"Warning: Only {nfft} points are being used for this spectra.")
    
    # Create a Hanning window of length nfft
    wind = np.hanning(nfft)
    # Normalization factor
    W1 = 2 / (np.linalg.norm(wind)**2)
    
    # First segment
    segment = detrend(x[:nfft])
    total = np.fft.fft(segment * wind)
    # Determine number of frequency points (MATLAB takes indices 2:floor(nfft/2+1))
    nfft_half = nfft // 2
    powe = total[1:nfft_half+1] * np.conjugate(total[1:nfft_half+1])
    
    # Loop over remaining segments if any (using fixed step size)
    if repeats > 1:
        step = int((max_ind - nfft) / (repeats - 1))
        for i in range(step, max_ind - nfft + 1, step):
            segment = detrend(x[i:i+nfft])
            total = np.fft.fft(segment * wind)
            powe += total[1:nfft_half+1] * np.conjugate(total[1:nfft_half+1])
    
    # Normalize PSD: note that in MATLAB the division is (W1 * powe) / (repeats * fs)
    powe = W1 * powe / (repeats * fs)
    
    # Create frequency vector: linearly spaced between fs/nfft and fs/2 with nfft/2 points.
    fre = np.linspace(fs/nfft, fs/2, nfft_half)
    
    return np.real(powe), fre

def stat_wave_1D(data, fs, fc, f_min, f_max, make_figure = True):
    """
    Process a 1D acceleration signal to compute displacement spectra, dominant frequency,
    significant wave height, and other wave statistics.
    
    Parameters:
        data (np.array): Input data (e.g., vertical acceleration).
        fs (float): Sampling frequency (Hz) of the input data.
        fc (float): Cutoff frequency (Hz) for high-pass filtering (to remove low frequencies).
        f_min (float): Minimum frequency (Hz) for significant wave height calculation.
        f_max (float): Maximum frequency (Hz) for significant wave height calculation.
        
    Returns:
        P_eta (np.array): Displacement Power Spectral Density (PSD).
        k (np.array): Frequency vector (Hz) corresponding to PSD values.
        a_z_filtered (np.array): High-pass filtered vertical acceleration.
        H_s (float): Significant wave height (meters).
        f_peak (float): Dominant (peak) frequency (Hz).
        T_peak (float): Peak period (seconds).
        t (np.array): Time vector (seconds) corresponding to the input data.
    """
    # 4th-order Butterworth high-pass filter
    order = 4
    Wn = fc / (fs / 2)  # Normalize cutoff frequency
    b, a = butter(order, Wn, btype='high')
    a_z_filtered = filtfilt(b, a, data)
    a_z = data  # Original acceleration data (unused later)
    t = np.arange(len(a_z)) / fs
    
    # Spectra calculation
    # Note: The MATLAB code creates a frequency vector "f" but does not use it further.
    # We now compute the PSD using the fast_psd function.
    # (MATLAB uses 2^11 for nfft and a PSD sampling frequency of 10 Hz.)
    P, k = fast_psd(a_z_filtered, 2**11, 10)
    
    # Convert acceleration PSD to displacement PSD:
    # P_eta = P / ((2*pi*k)^4), taking care to avoid division by zero.
    P_eta = P / ((2 * np.pi * k)**4)
    P_eta[k == 0] = 0  # Set the zero-frequency term to zero
    
    # Estimate dominant frequency and significant wave height
   
    idx = np.where((k >= f_min) & (k <= f_max))[0]
    if len(idx) == 0:
        print("No frequency components found in the specified range.")
        return None
    
    # Find the index of the maximum displacement PSD in the selected frequency range
    peak_idx_relative = np.argmax(P_eta[idx])
    peak_idx = idx[peak_idx_relative]
    f_peak = k[peak_idx]       # Dominant frequency (Hz)
    T_peak = 1 / f_peak        # Peak period (seconds)
    
    # Define a frequency range around the peak (±20% of f_peak)
    freq_range = [f_peak * 0.8, f_peak * 1.2]
    freq_indices = np.where((k >= freq_range[0]) & (k <= freq_range[1]))[0]
    
    # Integrate the displacement PSD over the frequency range
    m0 = np.trapz(P_eta[freq_indices], k[freq_indices])
    
    # Calculate Significant Wave Height (H_s)
    H_s = 4 * np.sqrt(m0)
    
    print(f"Significant Wave Height Hs: {H_s} meters")
    print(f"Dominant Frequency: {f_peak} Hz")
    print(f"Peak Period: {T_peak} seconds")
    

    if make_figure:
        # Plotting the PSD results
        plt.figure()
        plt.loglog(k, P, color=[0.2, 0.2, 0.2], linewidth=1.5, label='Acceleration PSD')
        plt.loglog(k, P_eta, 'b', linewidth=2, label='Displacement PSD')
        plt.plot(f_peak, P_eta[peak_idx], 'o', color='b', markersize=8, markerfacecolor='b')
        plt.axvline(x=freq_range[0], color='m', linestyle='--', linewidth=1.5)
        plt.axvline(x=freq_range[1], color='m', linestyle='--', linewidth=1.5)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power Spectral Density (m^2/Hz)')
        plt.legend()
        plt.show()
    
    return {
    'P_eta': P_eta,
    'k': k,
    'a_z_filtered': a_z_filtered,
    'H_s': H_s,
    'f_peak': f_peak,
    'T_peak': T_peak,
    't': t
}


#Example usage:
# data = np.loadtxt('your_data_file.txt')  # load your 1D data array
# fs = 10.0   # e.g., 100 Hz sampling frequency
# fc = 0.05     # e.g., 0.1 Hz cutoff frequency
# results = stat_wave_1D(data, fs, fc)
