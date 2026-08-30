function inst_freq = compute_instantaneous_frequency_analytical_fourier_series(amps, angular_freqs, t)
% COMPUTE_ANALYTICAL_INST_FREQ_FOURIER_SERIES analytically calculates the instantaneous frequency
% of a Fourier series E(t) = sum(amps_k * exp(1i * angular_freqs_k * t)) -
% which must not have +c.c
%
% INPUTS:
%   amps          : Vector of complex amplitudes
%   angular_freqs : Vector of angular frequencies
%   t             : Vector of time points
%
% OUTPUTS:
%   inst_freq     : Instantaneous frequency evaluated at time t

    % Ensure angular_freqs is a column vector
    angular_freqs = angular_freqs(:);
    
    % If amps is a vector, make it a column vector
    if isvector(amps)
        amps = amps(:);
    end
    
    % Ensure t is a row vector
    t_row = t(:).'; 
    
    % Compute the individual mode contributions: E_k = A_k * exp(1i * omega_k * t)
    % This uses implicit expansion if sizes are compatible
    E_components = amps .* exp(1i * angular_freqs .* t_row);
    
    % Sum over all modes to get the total field E
    E_tot = sum(E_components, 1);
    
    % The derivative of phase is Re[ sum( omega_k * E_k ) / E_tot ]
    numerator = sum(angular_freqs .* E_components, 1);
    
    inst_freq = real(numerator .* conj(E_tot)) ./ (abs(E_tot).^2 + eps);
    
    % Reshape output to match the shape of the input time vector t, only if amps was a single set of coefficients
    if isvector(t) && size(amps, 2) == 1
        inst_freq = reshape(inst_freq, size(t));
    end
end
