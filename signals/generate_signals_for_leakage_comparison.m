function [SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS] = generate_signals_for_leakage_comparison(freq_scaling, amp_scaling)

% --- 1. Base Parameters ---
angular_freqs_SO = [1, 2, 3, 4, 5] * 0.18 * freq_scaling;

amps_SO = [-0.156067704462866 + 0.331660754319902i,...
    -0.861836830340772 - 1.041781767817215i,...
    2.340666531884434 - 0.600981019561177i,...
    -0.502399901009243 + 2.633672512223463i,...
    -1.820362096071554 - 1.322570479164972i] * amp_scaling;

angular_freqs_COS = 1 * freq_scaling;

t_0 = 250;
T = 100;

% --- 2. Calculate In-Band Energy of the SO Signal ---
% By Parseval's theorem: Energy = (1/2pi) * integral( |X(w)|^2 ) dw

% Assuming the target bounds scale with the base frequency scaling:
w_lower = -1.05 * freq_scaling;
w_upper = -0.95 * freq_scaling;

% Format arrays into column vectors for MATLAB implicit expansion
w_n = angular_freqs_SO(:);
A_n = amps_SO(:);

% Define the Energy Spectral Density (ESD) function |X_SO(w)|^2
% The Fourier transform of a Gaussian-windowed exponential is a shifted Gaussian.
% This anonymous function is fully vectorized so `integral` can pass arrays of frequencies (w).
ESD_SO = @(w) abs( sum( conj(A_n) .* ...
    (T * sqrt(pi) * exp(-((w + w_n).^2) * T^2 / 4) .* exp(-1i * (w + w_n) * t_0)), 1 ) ).^2;

% Component Breakdown:
% abs(...).^2         : Energy Spectral Density definition (Magnitude Squared).
% sum(..., 1)         : Sums along the columns. MATLAB uses implicit expansion
%                       to combine the 5x1 (A_n, w_n) arrays with the 1xM (w)
%                       array, outputting a 5xM matrix. This sum collapses
%                       the 5 frequency components into one continuous spectrum.
% conj(A_n)           : Complex conjugate of the amplitude weights.
%
% --- The Analytical Fourier Transform Equation ---
% T * sqrt(pi)        : Base magnitude of the FT of a Gaussian envelope.
% exp(-((w+w_n).^2)...: The shape of the signal in the frequency domain.
%                       '(w + w_n)' represents the Frequency Shifting Property,
%                       moving the peak of the Gaussian to -w_n.
%                       'T^2 / 4' determines the variance/width of the spectrum.
% .* exp(-1i * ...)   : Phase shift caused by the Time Shifting Property. The
%                       envelope delay (t_0) introduces a linear phase rotation.
% =========================================================================

% Compute the in-band energy using numerical integration of the analytical FT
E_SO_in_band = (1 / (2 * pi)) * integral(ESD_SO, w_lower, w_upper);

% --- 3. Normalize the Cosine Signal ---
% We compute the total analytical time-domain energy of an UNNORMALIZED Cos signal.
% Unnormalized signal: x(t) = -exp(-i*w*t) * exp(-(t-t0)^2 / T^2)
% Int( |x(t)|^2 ) = Int( exp(-2(t-t0)^2 / T^2) ) = T * sqrt(pi / 2)
E_cos_unnormalized = T * sqrt(pi / 2);

% Since Energy is proportional to Amplitude^2, we scale using the square root.
% This sets the total energy of the Cos signal equal to the in-band energy of the SO signal.
Cos_normalization = sqrt(E_SO_in_band / E_cos_unnormalized);

% --- 4. Define Output Function Handles ---

% Function handle for cosine
Cos_signal = @(t) -conj(Cos_normalization .* exp(1i * angular_freqs_COS .* t) .* exp(-(t-t_0).^2./T^2));

% Function handle for superoscillatory signal (Expects 't' to be passed as a column vector)
SO_signal = @(t)  conj(sum( amps_SO .* exp(1i * angular_freqs_SO .* t), 2) .* exp(-(t-t_0).^2./T^2));

DO NUMERICAL VERIVICATION (UNIT TEST!!!!!) OF THE INTEGRAL USING FFT! 

end

% %% --- Script to Calculate Signal Energy Ratio ---
% 
% % 1. Set scaling parameters (you can adjust these as needed)
% freq_scaling = 1;
% amp_scaling = 1;
% 
% % 2. Retrieve the signal function handles
% [SO_signal, Cos_signal, ~, ~] = generate_signals_for_leakage_comparison(freq_scaling, amp_scaling);
% 
% % 3. Define the Power wrappers (Magnitude Squared)
% % We transpose 't' (using .') because SO_signal expects a column vector,
% % and we transpose the result back so 'integral' gets the row vector it expects.
% power_SO  = @(t) (abs(SO_signal(t.')).^2).';
% power_Cos = @(t) (abs(Cos_signal(t.')).^2).';
% 
% % 4. Compute the total energy by integrating from -inf to inf
% % (The Gaussian window exp(-(t-t_0)^2/T^2) ensures this converges quickly)
% energy_SO  = integral(power_SO, -inf, inf);
% energy_Cos = integral(power_Cos, -inf, inf);
% 
% % 5. Calculate the ratio
% energy_ratio = energy_Cos/ energy_SO;
% 
% % 6. Display the results
% fprintf('--- Energy Comparison ---\n');
% fprintf('Energy of SO signal:  %.4f\n', energy_SO);
% fprintf('Energy of Cos signal: %.4f\n', energy_Cos);
% fprintf('Ratio (Cos / So):     %.4e\n', energy_ratio);