function [] = plot_signals(VinSOFun, VinCosFun, angular_freqs,angular_freqs_COS)

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

figure
hold on;
% Capture the handles (h1, h2) as you plot
h2 = plot(t_axis, real(sampled_signal), '-','color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'DisplayName', '\textbf{SO}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');
legend([h1, h2]);

PlotUtils.styleAxes(gca);
hold off;


%% FFT the signals
% Compute FFT
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

fft_cos = fftshift(abs(fft(sampled_signal,N))) / N;
fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N))) / N;


%% Plot FFTs
figure;
hold on;

% 1. Plotting Data
h2 = plot(freq_axis, (fft_cos), '-', 'Color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(freq_axis, fft_superoscillation, '-', 'Color', 'r', 'DisplayName', '\textbf{SO}');

xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');
legend([h1, h2]);

PlotUtils.styleAxes(gca);
hold off;

end