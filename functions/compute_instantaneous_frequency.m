function [inst_freq, phase, z] = compute_instantaneous_frequency(y, t_or_dt)
%COMPUTE_INSTANTANEOUS_FREQUENCY Analytically calculates the instantaneous frequency of a signal.
%
%   [inst_freq, phase, z] = compute_instantaneous_frequency(y, t_or_dt)
%
%   INPUTS:
%     y       : Vector of signal values (real or complex). If complex, the real part is used.
%     t_or_dt : Either a vector of time points corresponding to y, or a scalar sampling interval (dt).
%
%   OUTPUTS:
%     inst_freq : Instantaneous frequency of the signal.
%     phase     : Unwrapped phase of the analytic signal.
%     z         : Complex analytic signal.


    % Get the real part of the signal as the Hilbert transform is defined on real signals
    y_real = real(y);

    % Calculate the analytic signal
    z = hilbert(y_real);

    % Calculate the unwrapped phase
    phase = unwrap(angle(z));

    % Compute the gradient (derivative) of phase with respect to time/sampling step
    inst_freq = gradient(phase, t_or_dt);

end
