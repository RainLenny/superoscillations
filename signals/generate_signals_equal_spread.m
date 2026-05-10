function [SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS] = generate_signals_equal_spread(freq_scaling, amp_scaling)

% Parameters
angular_freqs_SO = [1,2,3,4,5] * 0.18 * freq_scaling;

amps_SO    = [-0.156067704462866 + 0.331660754319902i,...
    -0.861836830340772 - 1.041781767817215i,...
    2.340666531884434 - 0.600981019561177i,...
    -0.502399901009243 + 2.633672512223463i,...
    -1.820362096071554 - 1.322570479164972i] * amp_scaling;

% Contain the same frequency components
angular_freqs_COS = angular_freqs_SO;

t_0 = 250;
T = 100;

% Function handles
SO_signal = @(t) reshape(conj(sum( amps_SO .* exp(1i * angular_freqs_SO .* t(:)), 2) .* exp(-(t(:)-t_0).^2./T^2)), size(t));

amps_COS_unscaled = [1,2,3,4,5];
Cos_signal_unscaled = @(t) reshape(conj(sum( amps_COS_unscaled .* exp(1i * angular_freqs_COS .* t(:)), 2) .* exp(-(t(:)-t_0).^2./T^2)), size(t));

% Integration bounds
int_min = -inf;
int_max = inf;

% Calculate energies
energy_SO = integral(@(t) abs(SO_signal(t)).^2, int_min, int_max);
energy_Cos_unscaled = integral(@(t) abs(Cos_signal_unscaled(t)).^2, int_min, int_max);

% Apply scaling factor
scale_factor = sqrt(energy_SO / energy_Cos_unscaled);
amps_COS = amps_COS_unscaled * scale_factor;

% Final Cosine Signal
SO_signal = @(t) reshape(conj(sum( amps_SO   .* exp(1i * angular_freqs_SO  .* t(:)), 2) .* exp(-(t(:)-t_0).^2./T^2)), size(t));

Cos_signal = @(t) reshape(conj(sum( amps_COS .* exp(1i * angular_freqs_COS .* t(:)), 2) .* exp(-(t(:)-t_0).^2./T^2)), size(t));

%% INBUILT VERIFICATION (UNIT TEST)

% Calculate the actual energy of the newly scaled Cos_signal
energy_Cos_final = integral(@(t) abs(Cos_signal(t)).^2, int_min, int_max);

% Check if they match within a tight numerical tolerance
tolerance = 1e-8;
if abs(energy_SO - energy_Cos_final) > tolerance
    % Throw a custom error breaking execution if the math failed
    error('SignalGen:EnergyMismatch', ...
        'Validation Failed: Energies do not match.\nEnergy SO: %.10f\nEnergy COS: %.10f\nDifference: %.10f', ...
        energy_SO, energy_Cos_final, abs(energy_SO - energy_Cos_final));
end

end