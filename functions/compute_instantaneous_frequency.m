function inst_freq = compute_instantaneous_frequency(angular_freqs, amps, t)
%COMPUTE_INSTANTANEOUS_FREQUENCY Analytically calculates the instantaneous frequency of a signal
%consisting of a sum of complex exponentials.
%
%   inst_freq = compute_instantaneous_frequency(angular_freqs, amps, t)
%
%   INPUTS:
%     angular_freqs : Vector of angular frequencies
%     amps          : Vector of complex amplitudes
%     t             : Vector of time points
%
%   OUTPUTS:
%     inst_freq : Instantaneous frequency of the signal evaluated at points t.

    angular_freqs = angular_freqs(:);
    amps = amps(:);
    t_shape = size(t);
    t = t(:).'; % make row vector
    
    % Denominator: sum(c_n * exp(i * w_n * t))
    f_t = sum(amps .* exp(1i * angular_freqs .* t), 1);
    
    % Numerator: sum(c_n * w_n * exp(i * w_n * t))
    df_t = sum(amps .* angular_freqs .* exp(1i * angular_freqs .* t), 1);
    
    % Instantaneous frequency is Re(df_t / f_t)
    inst_freq = real(df_t ./ f_t);
    
    % Return same shape as input t
    inst_freq = reshape(inst_freq, t_shape);

end
