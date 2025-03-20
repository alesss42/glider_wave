function J = costfunction(D, S_meas_k, S1D, RAO_heave, RAO_pitch, RAO_roll, DeltaTheta)
    % D is a column vector of length N_dir representing D(theta) at this frequency.
    % S_meas_k is the measured 3x3 cross-spectral matrix at this frequency.
    % S1D is the scalar 1D wave spectrum at this frequency.
    % RAO_heave, RAO_pitch, RAO_roll are row vectors (1 x N_dir) of RAO values.
    % DeltaTheta is the angular bin width.
    
    % Compute theoretical auto-spectra and cross-spectra
    S11 = S1D * sum(RAO_heave.^2 .* D') * DeltaTheta;
    S22 = S1D * sum(RAO_pitch.^2 .* D') * DeltaTheta;
    S33 = S1D * sum(RAO_roll.^2  .* D') * DeltaTheta;
    S12 = S1D * sum(RAO_heave .* RAO_pitch .* D') * DeltaTheta;
    S13 = S1D * sum(RAO_heave .* RAO_roll  .* D') * DeltaTheta;
    S23 = S1D * sum(RAO_pitch .* RAO_roll  .* D') * DeltaTheta;
    
    % Construct the theoretical cross-spectral matrix
    S_theo = [S11, S12, S13; 
              conj(S12), S22, S23; 
              conj(S13), conj(S23), S33];
          
    % Calculate cost as the squared Frobenius norm of the difference
    diffMat = S_meas_k - S_theo;
    J = norm(diffMat, 'fro')^2;
end
