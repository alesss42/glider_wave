function [P_eta, k,a_z_filtered, H_s, f_peak, T_peak, t] = stat_wave_1D(data, fs, fc)
    % data: input data (e.g., glider_data.ACC_0_)
    % fs: sampling frequency in Hz
    % fc: cutoff frequency in Hz (removing low frequencies)


    % Data preparation 
    % 4th-order Butterworth high-pass filter
    order = 4;  
    Wn = fc / (fs/2); % Normalize frequency
    [b, a] = butter(order, Wn, 'high');

    a_z_filtered = filtfilt(b, a, data); % Filtered vertical acceleration data
    a_z = data; % Original vertical acceleration data
    t = (0:length(a_z)-1) / fs; % Time vector

    


    % Spectra calculation 
    N = length(a_z_filtered); 
    dt = 1/fs; % Time step
    f = (0:N/2-1) * (fs/N); 

    % Computes the Power Spectral Density (PSD)
    [P, k] = fast_psd(a_z_filtered, 2^11, 10); % 2^10 is the number of points in the FFT
    
    % loglog(k,P); hold on
    % Convert acceleration PSD to displacement PSD
    P_eta = P ./ ( (2 * pi * k).^4 ); 
    P_eta(k == 0) = 0; 

    
    % loglog(k, P_eta, 'b', 'LineWidth', 2);

    % Estimate dominant frequency and significant wave height

    % Here we need to figure out better what our boudnaries would look like,
    f_min = 0.05; % Minimum frequency to consider
    f_max = 10.0; % Maximum frequency to consider
    idx = find(k >= f_min & k <= f_max);
    [~, peak_idx] = max(P_eta(idx)); % Find max only in this range

    % Find the peak frequency and period
    f_peak = k(idx(peak_idx)); % Peak frequency
    T_peak = 1 / f_peak;  % Peak period

    % Frequency range around the peak (±20% of the peak frequency)
    freq_range = [f_peak * 0.8, f_peak * 1.2];  % Range around the peak frequency
    freq_indices = find(k >= freq_range(1) & k <= freq_range(2));

    % % Integrate over the selected frequency range (S_eta over the range)
     m0 = trapz(k(freq_indices), P_eta(freq_indices));  % Numerical integration over the selected range

    % % Calculate Significant Wave Height (H_s)
     H_s = 4 * sqrt(m0);

    % % Display result
     disp(['Significant Wave Height Hs: ', num2str(H_s), ' meters']);
     disp(['Dominant Frequency: ', num2str(f_peak), ' Hz']);
     disp(['Peak Period: ', num2str(T_peak), ' seconds']);

    % % Figure

    loglog(k, P, 'k', 'color', [.2 .2 .2] ,'LineWidth', 1.5); hold on;
    loglog(k, P_eta, 'b', 'LineWidth', 2);
    line(f_peak, P_eta(idx(peak_idx)), 'Marker', 'o', 'Color', 'b', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    xline(freq_range(1), 'm--', 'LineWidth', 1.5);
    xline(freq_range(2), 'm--', 'LineWidth', 1.5);

    ylabel('Power Spectral Density (m^2/Hz)');

end
