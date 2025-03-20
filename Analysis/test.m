% Read the file and skip the first row explicitly
opts = detectImportOptions('../data/Regular_wave_GLIDER.csv');
opts.DataLines = [2 Inf]; % Start reading from the second row to the end
data = readtable('../data/Regular_wave_GLIDER.csv', opts);


opts = detectImportOptions('../data/JONSWAP_GLIDER.csv');
opts.DataLines = [2 Inf]; % Start reading from the second row to the end
data2 = readtable('../data/JONSWAP_GLIDER.csv', opts);

%% figure
figure
[P,k]=fast_psd(data2.X,2^10,1/nanmean(diff(data2.Time)));
loglog(k,P );

% 
%shg; hold on; text(.25,1e-5,‘Surface waves’); text(.04,.9e-6,‘-5/3’,‘color’,‘r’); text(2e-3,1,‘-3’,‘color’,‘b’); hold off

hold on
[P,k]=fast_psd(data2.Y,2^10,1/nanmean(diff(data2.Time)));
loglog(k,P );

hold on
[P,k]=fast_psd(data2.Z,2^10,1/nanmean(diff(data2.Time)));
loglog(k,P );
legend('X', 'Y', 'Z')
xlabel('Frequency')
ylabel('Power')


%%
figure
[P,k]=fast_psd(data.X,2^8,1/nanmean(diff(data.Time)));
loglog(k,P );
% 
%shg; hold on; text(.25,1e-5,‘Surface waves’); text(.04,.9e-6,‘-5/3’,‘color’,‘r’); text(2e-3,1,‘-3’,‘color’,‘b’); hold off

hold on
[P,k]=fast_psd(data.Y,2^8,1/nanmean(diff(data.Time)));
loglog(k,P );

hold on
[P,k]=fast_psd(data.Z,2^8,1/nanmean(diff(data.Time)));
loglog(k,P );

%%
fs = 100; % Sampling frequency (Hz)
t = 0:1/fs:10-1/fs; % Time vector (10 seconds)
y = 3*sin(2*pi*2*t) + 2*cos(2*pi*5*t) + 0.5*randn(size(t)); % Signal with noise

frequencies = [2]; % Known frequencies (Hz)

X = []; % Design matrix
for f = frequencies
    X = [X, cos(2*pi*f*t)', sin(2*pi*f*t)'];
end

coeffs = X \ y'; % Solve using least squares
amplitudes = sqrt(coeffs(1:2:end).^2 + coeffs(2:2:end).^2);
phases = atan2(-coeffs(2:2:end), coeffs(1:2:end));

y_reconstructed = X * coeffs;

figure;
plot(t, y, 'k', 'DisplayName', 'Original Signal'); hold on;
plot(t, y_reconstructed, 'r--', 'DisplayName', 'Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend;
title('Time Series Decomposition into Known Frequencies');
grid on;
