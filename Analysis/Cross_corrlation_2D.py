import numpy as np
from scipy.signal import butter, filtfilt
from scipy.signal import csd, coherence

# Define the window and overlap
window = np.hanning(Nfft)
noverlap = int(Nfft / 5)  # 20% overlap


def quick_high_pass_filter(data, fs, fc, order):
    """
    Apply a high-pass Butterworth filter to the input data.

    Parameters:
        data (array_like): Input signal data.
        fs (float): Sampling frequency in Hz.
        fc (float): Cutoff frequency in Hz.
        order (int): Order of the Butterworth filter.

    Returns:
        a_z_filtered (ndarray): The high-pass filtered signal.
    """
    # Normalize the cutoff frequency
    Wn = fc / (fs / 2)
    
    # Design the Butterworth high-pass filter
    b, a = butter(order, Wn, btype='high')
    
    # Apply zero-phase filtering using filtfilt
    a_z_filtered = filtfilt(b, a, data)
    
    return a_z_filtered


# Example: glider_data is assumed to be a dictionary with NumPy arrays.
# For example:
# glider_data = {
#     'ACC_0_': np.array([...]),
#     'ACC_1_': np.array([...]),
#     'ACC_2_': np.array([...]),
#     'PQR_0_': np.array([...]),
#     'PQR_1_': np.array([...]),
#     'PQR_2_': np.array([...]),
#     'Euler_0_': np.array([...]),
#     'Euler_1_': np.array([...]),
#     'Euler_2_': np.array([...])
# }
# And thr is defined as some threshold value.

# Preparing the data: create a boolean mask where the condition is met.
ac = glider_data['ACC_0_'] < thr

# Create an empty dictionary to store the filtered glider data.
glider = {}

# Acceleration (m/s^2)
glider['heave'] = glider_data['ACC_0_'][ac]   # Heave (x) up-down motion
glider['surge'] = glider_data['ACC_1_'][ac]    # Surge (y) forward-backward motion
glider['sway']  = glider_data['ACC_2_'][ac]    # Sway (z) side-to-side motion

# Angular velocity (rad/s)
glider['PQR_0'] = glider_data['PQR_0_'][ac]
glider['PQR_1'] = glider_data['PQR_1_'][ac]
glider['PQR_2'] = glider_data['PQR_2_'][ac]

# Euler angles (deg)
glider['yaw']   = glider_data['Euler_0_'][ac]   # Yaw
glider['pitch'] = glider_data['Euler_1_'][ac]   # Pitch
glider['roll']  = glider_data['Euler_2_'][ac]   # Roll

# Now, glider contains the filtered data.
glider['pitch_rad'] = np.deg2rad(glider['pitch'])
glider['roll_rad']  = np.deg2rad(glider['roll'])


import numpy as np

# Set the number of FFT points
fs = 10  # Sampling frequency (Hz)
fc = 0.05  # Cutoff frequency (Hz)
order = 4  # Filter order
# Define the number of FFT points (Nfft)
Nfft = 2**10  # Number of FFT points

# Choose a window size and overlap for mscohere (if needed in your implementation)
window = np.hanning(Nfft)
noverlap = Nfft // 2  # 50% overlap

# Prepare dictionaries to hold the filtered data and PSD results
PSD = {}
data_p = {}

# Iterate over all keys (fields) in the glider dictionary
for var_name, data in glider.items():
    # High-pass filter: quick_high_pass_filter should be defined elsewhere
    data_p[var_name] = quick_high_pass_filter(data, fs, fc, order)  # Removing low frequency noise
    
    # Compute the Power Spectral Density (PSD) for the filtered data
    PSD[var_name], k_psd = fast_psd(data_p[var_name], Nfft, fs)  # Similar to MATLAB's pwelch



# Compute Cross-Spectral Density (CSD)

# Heave-Pitch CSD
f_csd, CSD_hp = csd(data_p['heave'], data_p['pitch_rad'], fs=fs, window=window, noverlap=noverlap, nfft=Nfft)
k_csd = f_csd  # Frequency vector for CSD

# Heave-Roll CSD
_, CSD_hr = csd(data_p['heave'], data_p['roll_rad'], fs=fs, window=window, noverlap=noverlap, nfft=Nfft)

# Pitch-Roll CSD
_, CSD_pr = csd(data_p['pitch'], data_p['roll_rad'], fs=fs, window=window, noverlap=noverlap, nfft=Nfft)

# Compute Magnitude-Squared Coherence

# Heave-Pitch Coherence
f_coh, Coherence_hp = coherence(data_p['heave'], data_p['pitch_rad'], fs=fs, window=window, noverlap=noverlap, nfft=Nfft)

# Heave-Roll Coherence
_, Coherence_hr = coherence(data_p['heave'], data_p['roll_rad'], fs=fs,  window=window, noverlap=noverlap, nfft=Nfft)

