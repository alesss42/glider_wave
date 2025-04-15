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