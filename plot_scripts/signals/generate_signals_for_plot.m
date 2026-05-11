function [SO_signal, Cos_signal,angular_freqs_SO, angular_freqs_COS] = generate_signals_for_plot(freq_scaling,amp_scaling)

% Parameters
load('SO_signal_data.mat', 'amps_SO', 'angular_freqs_SO')

angular_freqs_COS = 0.9 * freq_scaling; %angular_freqs_SO(end);

t_0 = 250;
T = 100;

Cos_normalization = 4.4781551283; %Scaling such that the energy of the actual drives would be the same

% Function handle for cosine
Cos_signal = @(t)   -real(Cos_normalization .* cos(angular_freqs_COS .* t).*exp(-(t-t_0).^2./T^2));

% Function handle for superoscillatory signal
SO_signal = @(t)  real(sum( amps_SO .* exp(1i * angular_freqs_SO .* t),2).*exp(-(t-t_0).^2./T^2));
end