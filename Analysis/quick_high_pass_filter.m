function a_z_filtered = quick_high_pass_filter(data, fs, fc, order)

    % Design a 4th-order Butterworth high-pass filter
Wn = fc / (fs/2); % Normalize frequency
[b, a] = butter(order, Wn, 'high');

% Apply the filter to the vertical acceleration data
a_z_filtered = filtfilt(b, a, data); % Assuming glider_data.Var9 is your vertical acceleration data
end