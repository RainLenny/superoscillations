function [SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS] = generate_signals_Baranov(freq_scaling, amp_scaling)

% Parameters
angular_freqs_SO = [1,2,3,4,5] * 0.18 * freq_scaling;

amps_SO    = [-0.156067704462866 + 0.331660754319902i,...
    -0.861836830340772 - 1.041781767817215i,...
    2.340666531884434 - 0.600981019561177i,...
    -0.502399901009243 + 2.633672512223463i,...
    -1.820362096071554 - 1.322570479164972i] * amp_scaling;

angular_freqs_COS = 0.9 * freq_scaling;

t_0 = 250;
T = 100;

% Function handle for superoscillatory signal
% (Note: works best when t is a column vector due to sum(..., 2))
SO_signal = @(t) conj(sum( amps_SO .* exp(1i * angular_freqs_SO .* t), 2) .* exp(-(t-t_0).^2./T^2));

% --- Energy Calculation ---
% MATLAB's 'integral' passes t as a row vector. We use a wrapper t(:) inside
% the function to force it into a column vector to satisfy your SO_signal's
% matrix dimensions, and then reshape it back to a row vector for 'integral'.

SO_power_integrand = @(t) reshape(abs(SO_signal(t(:))).^2, size(t));
E_SO = integral(SO_power_integrand, -inf, inf);

% Calculate the energy of the unnormalized cosine envelope
Cos_unnorm_integrand = @(t) reshape(abs(exp(-(t(:)-t_0).^2./T^2)).^2, size(t));
E_Cos_unnorm = integral(Cos_unnorm_integrand, -inf, inf);

% Dynamically calculate normalization to match the energy of SO_signal
% (Energy scales with the square of the amplitude, so we take the square root)
Cos_normalization = sqrt(E_SO / E_Cos_unnorm);

% Function handle for cosine using the dynamic normalization
Cos_signal = @(t) -conj(Cos_normalization .* exp(1i * angular_freqs_COS .* t) .* exp(-(t-t_0).^2./T^2));

%% Numerical Test
% Calculate the energy of the newly normalized Cosine signal
Cos_power_integrand = @(t) reshape(abs(Cos_signal(t(:))).^2, size(t));
E_Cos = integral(Cos_power_integrand, -inf, inf);

% Verify the energies are the same within a small numerical tolerance
tol = 1e-6;
if abs(E_SO - E_Cos) > tol
    error('Energy normalization failed! Difference exceeds tolerance.\n E_SO: %g\n E_Cos: %g\n Diff: %g', ...
        E_SO, E_Cos, abs(E_SO - E_Cos));
end

end