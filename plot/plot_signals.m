function [] = plot_signals(VinSOFun, VinCosFun, angular_freqs, angular_freqs_COS, Flat_signal, Rand_signal)

% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();

%% SAMPLING Constants
f_sampling = 60/(2*pi);
% fundamental_period of the superoscillating signal:
T_period = compute_fundamental_period([angular_freqs,angular_freqs_COS],f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample‐interval in seconds
t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1);
t_axis = t_axis(:);


%% Sample the signals
sampled_signal = VinCosFun(t_axis);
sampled_superoscillation = VinSOFun(t_axis);

has_four = (nargin >= 6 && ~isempty(Flat_signal) && ~isempty(Rand_signal));
if has_four
    sampled_flat = Flat_signal(t_axis);
    sampled_rand = Rand_signal(t_axis);
end

figure
hold on;
% Capture the handles as you plot
h2 = plot(t_axis, real(sampled_signal), '-','color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'DisplayName', '\textbf{SO}');
if has_four
    h3 = plot(t_axis, real(sampled_flat), '-','color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');
    h4 = plot(t_axis, real(sampled_rand), '-','color', 'm', 'DisplayName', '\textbf{Rand Phase}');
    legend([h1, h2, h3, h4]);
else
    legend([h1, h2]);
end

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');

PlotUtils.styleAxes(gca);
hold off;


%% FFT the signals
% Compute FFT
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

fft_cos = fftshift(abs(fft(sampled_signal,N))) / N;
fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N))) / N;
if has_four
    fft_flat = fftshift(abs(fft(sampled_flat,N))) / N;
    fft_rand = fftshift(abs(fft(sampled_rand,N))) / N;
end


%% Plot FFTs
figure;
hold on;

% 1. Plotting Data
h2 = plot(freq_axis, (fft_cos), '-', 'Color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(freq_axis, fft_superoscillation, '-', 'Color', 'r', 'DisplayName', '\textbf{SO}');
if has_four
    h3 = plot(freq_axis, fft_flat, '-', 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');
    h4 = plot(freq_axis, fft_rand, '-', 'Color', 'm', 'DisplayName', '\textbf{Rand Phase}');
    legend([h1, h2, h3, h4]);
else
    legend([h1, h2]);
end

xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');

PlotUtils.styleAxes(gca);
hold off;

end