function [SO_signal, Cos_signal,angular_freqs_SO, angular_freqs_COS] = generate_signals_for_plot(freq_scaling,amp_scaling)

% Parameters
angular_freqs_SO = [1,2,3,4,5] * 0.18 * freq_scaling;

amps_SO    = [-0.156067704462866 + 0.331660754319902i,...
    -0.861836830340772 - 1.041781767817215i,...
    2.340666531884434 - 0.600981019561177i,...
    -0.502399901009243 + 2.633672512223463i,...
    -1.820362096071554 - 1.322570479164972i] * amp_scaling;

angular_freqs_COS = 0.9 * freq_scaling; %angular_freqs_SO(end);

t_0 = 250;
T = 100;

Cos_normalization = 4.478155155550048; %Scaling such that the energy of the actual drives would be the same


% Function handle for cosine
Cos_signal = @(t)   real(Cos_normalization .* exp(1i * angular_freqs_COS .* t));

% Function handle for superoscillatory signal
SO_signal = @(t)  real(sum( amps_SO .* exp(1i * angular_freqs_SO .* t),2));
end