clear; clc;
%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();
[VinSOFun_comp, angular_freqs, amps_SO] = generate_SO_from_dat_file('2p5_cos_Derek', 1, 1, false);
[VinCosFun_comp, angular_freqs_COS] = generate_Cos_reference(1, 1, false);

VinSOFun = @(t) real(VinSOFun_comp(t)) ./ real(VinSOFun_comp(0));

VinCosFun = @(t) real(VinCosFun_comp(t));


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
% Using a cos reference of frequency 1 (\omega_0) as requested
cos1 = VinCosFun(t_axis);
sampled_superoscillation = VinSOFun(t_axis);


%% FFT the signals
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

% Calculate the FFT magnitude
fft_superoscillation = fftshift(abs(fft(sampled_superoscillation, N)))*dt;
fft_cos = fftshift(abs(fft(cos1, N)))*dt;

% Normalize by the number of samples (N) to get true Fourier coefficients
fft_superoscillation = fft_superoscillation / (sum(fft_superoscillation));
fft_cos = fft_cos / (sum(fft_cos));


%% Figure (a): The signals and the ideal filter in the time domain.
figure;
hold on;
% Scale filter to have similar amplitude to signals for visualization
filter_time = (0.7/pi) * sinc(0.7 * t_axis / pi);
filter_time = filter_time / max(filter_time) ; % scale factor to match image roughly

plot(t_axis, real(sampled_superoscillation), '-', 'color', 'r', 'DisplayName', '\textbf{SO}');
plot(t_axis, cos1, ':', 'color', 'b', 'DisplayName', '\textbf{\boldmath$\mathbf{\omega_0}$}');
plot(t_axis, filter_time, '-.', 'color', 'k', 'DisplayName', '\textbf{Filter}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');
legend('show');
xlim([-8, 8]);
PlotUtils.styleAxes(gca);
hold off;


%% Figure (b): The signals and the ideal filter on the frequency domain.
figure;
hold on;
% Filter in frequency domain: Rect function from -0.7 to 0.7
filter_freq = double(abs(freq_axis) <= 0.7);

% Assign handles to each plot
h1 = plot(freq_axis, fft_cos, '-o', 'Color', 'b', 'DisplayName', '\textbf{\boldmath$\mathbf{\omega_0}$}');
h2 = plot(freq_axis, fft_superoscillation, '-o', 'Color', 'r', 'DisplayName', '\textbf{SO}');
h3 = plot(freq_axis, filter_freq, '--', 'Color', 'k', 'DisplayName', '\textbf{Filter}');

xlabel('\boldmath$\mathbf{Frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [norm.]}$');

% Specify the legend order explicitly using the handles
legend([h2, h1, h3]);

xlim([0, 1.5]);
PlotUtils.styleAxes(gca);
hold off;


%% Figure (c): The signals filtered by the ideal filter.
% Analytical filtering:
% SO is untouched (bandwidth < 0.7)
filtered_superoscillation = sampled_superoscillation;
% Cos is completely filtered out (freq 1 > 0.7)
filtered_cos = zeros(size(cos1));

figure;
tlo = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Superoscillating Signal
ax1 = nexttile;
hold on;
plot(t_axis, real(sampled_superoscillation), '-', 'color', 'r', 'DisplayName', '\textbf{Original}');
plot(t_axis, real(filtered_superoscillation), '--', 'color', 'b', 'DisplayName', '\textbf{Filtered}');
title('\textbf{SO Signal}');
xlim([-12, 12]);
PlotUtils.styleAxes(gca);
hold off;

% Subplot 2: Cos Signal
ax2 = nexttile;
hold on;
plot(t_axis, cos1, '-', 'color', 'r', 'DisplayName', '\textbf{Original}');
plot(t_axis, filtered_cos, '--', 'color', 'b', 'DisplayName', '\textbf{Filtered}');
xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
title('\textbf{\boldmath$\mathbf{COS(\omega_0)}$}');
xlim([-12, 12]);
ylim([-1.2, 1.2]); % Adjust ylim for visibility of the zero line
PlotUtils.styleAxes(gca);
hold off;

lgd = legend(ax1, 'show', 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'south';

% Add shared y-axis label using native tiledlayout feature
ylab = ylabel(tlo, '\boldmath$\mathbf{Amplitude \ [arb.]}$');
ylab.FontSize = ax1.FontSize;
ylab.Interpreter = 'latex';
