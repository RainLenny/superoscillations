function T_period = compute_fundamental_period(angular_freqs,f_sampling)
% Compute the fundamental period (T_period) of a superoscillating signal
% Using its angular frequencies.

% Set tolerance for rational approximation (defines how close the rational
% approximation must be to the true value.)
tol = 1e-12;

% Convert angular frequencies (ω, in rad/s) to cycles per sample
f_cycles = angular_freqs / (2*pi * f_sampling);

% Find rational approximations for each normalized frequency
% rat(x, tol) returns numerator/denominator such that numerator/denominator ≈ x
% we only need the denominators to determine when each component repeats
[~, denoms] = rat(f_cycles, tol);

% Compute Least Common Multiple (LCM) of all denominators to get the
% fundamental period in samples
prediod_samples = denoms(1);
for d = denoms(2:end)
    prediod_samples = lcm(prediod_samples, d);
end

% Convert the period from samples back to seconds
% Divide by the sampling rate to get the actual time duration
T_period = prediod_samples / f_sampling;

end


