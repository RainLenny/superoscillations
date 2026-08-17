clear; clc;
%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();

[VinSOFun_comp, angular_freqs, amps_SO] = generate_SO_from_dat_file('2p5_cos_Derek', 1, 1, false);
VinSOFun = @(t) real(VinSOFun_comp(t));

%% Parameters
omega_0 = 1;
omega_c = 0.7; % Filter cutoff frequency

%% Time Axis Setup
f_sampling = 600/(2*pi);
dt = 1 / f_sampling;
% Use a fixed range for better visualization matching the image
t_axis = -15 : dt : 15; 
t_axis = t_axis(:);

%% Generate Time Domain Signals
so_signal = VinSOFun(t_axis);
% Normalize SO signal so peak is 1 (to match the visualization)
max_so = max(abs(so_signal));
if max_so > 0
    so_signal = so_signal / max_so;
end

% Cosine reference signal
cos_signal = cos(omega_0 * t_axis);

% Ideal filter impulse response in time domain: h(t) = (omega_c / pi) * sinc(omega_c * t / pi)
% MATLAB's sinc is sinc(x) = sin(pi*x)/(pi*x), so sinc(omega_c * t / pi) gives sin(omega_c * t) / (omega_c * t)
filter_time = (omega_c / pi) * sinc(omega_c * t_axis / pi);

%% Figure 1: Time Domain
figure;
hold on;
plot(t_axis, so_signal, '-', 'Color', '#FF6666', 'LineWidth', 4, 'DisplayName', '\textbf{Superoscillation}');
plot(t_axis, cos_signal, '--', 'Color', '#8888FF', 'LineWidth', 4, 'DisplayName', '\textbf{Cosine}');
plot(t_axis, filter_time, '-.', 'Color', [0.7 0.7 0.7], 'LineWidth', 4, 'DisplayName', '\textbf{Filter}');
hold off;

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb. \ units]}$');
title('\boldmath$\mathbf{(a) \ The \ signals \ and \ the \ ideal \ filter \ in \ the \ time \ domain.}$');
legend('show', 'Location', 'best');
xlim([-15, 15]);
ylim([-0.5, 1.1]);
set(gca, 'XColor', 'k');
PlotUtils.styleAxes(gca);

%% Figure 2: Frequency Domain
figure;
hold on;

% Filter frequency response (ideal rect from -omega_c to omega_c)
plot([-1.5, -omega_c, -omega_c, omega_c, omega_c, 1.5], ...
     [0, 0, 1, 1, 0, 0], ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 3, 'DisplayName', '\textbf{Filter}');

% SO frequencies (we know angular_freqs are < omega_c)
if ~isempty(angular_freqs)
    stem_freqs = [angular_freqs(:); -angular_freqs(:)];
    stem_amps = [abs(amps_SO(:)); abs(amps_SO(:))]; 
    % Normalize stems to 1 for visual clarity
    if max(stem_amps) > 0
        stem_amps = stem_amps / max(stem_amps);
    end
    stem(stem_freqs, stem_amps, 'filled', 'Color', '#FF6666', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '\textbf{Superoscillation}');
else
    % Fallback if empty
    stem([0], [1], 'filled', 'Color', '#FF6666', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '\textbf{Superoscillation}');
end

% Cosine frequency (deltas at -omega_0 and omega_0)
plot([-omega_0, -omega_0], [0, 1], '-', 'Color', '#8888FF', 'LineWidth', 3, 'DisplayName', '\textbf{Cosine}');
plot([omega_0, omega_0], [0, 1], '-', 'Color', '#8888FF', 'LineWidth', 3, 'HandleVisibility', 'off');

hold off;
xlabel('\boldmath$\mathbf{Frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb. \ units]}$');
title('\boldmath$\mathbf{(b) \ The \ signals \ and \ the \ ideal \ filter \ on \ the \ frequency \ domain.}$');
legend('show', 'Location', 'best');
xlim([-1.5, 1.5]);
ylim([0, 1.2]);
set(gca, 'XColor', 'k');
PlotUtils.styleAxes(gca);

%% Figure 3: Filtered Signals (Analytical)
figure;

% Analytical filtering:
% 1. SO is comprised of frequencies < omega_c, so it is untouched.
filtered_so = so_signal; 
% 2. Cosine has frequency omega_0 > omega_c, so it is completely filtered out.
filtered_cos = zeros(size(t_axis));

% Subplot 1: Superoscillating Signal
subplot(2,1,1);
hold on;
plot(t_axis, so_signal, '-', 'Color', '#FF6666', 'LineWidth', 4, 'DisplayName', '\textbf{Original}');
plot(t_axis, filtered_so, '--', 'Color', '#8888FF', 'LineWidth', 4, 'DisplayName', '\textbf{Filtered}');
hold off;
title('\boldmath$\mathbf{Superoscillating \ Signal}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb. \ units]}$');
legend('show', 'Location', 'best');
xlim([-15, 15]);
ylim([-0.5, 1.1]);
set(gca, 'XColor', 'k');
PlotUtils.styleAxes(gca);

% Subplot 2: Cosine Signal
subplot(2,1,2);
hold on;
plot(t_axis, cos_signal, '-', 'Color', '#FF6666', 'LineWidth', 4, 'DisplayName', '\textbf{Original}');
plot(t_axis, filtered_cos, '--', 'Color', '#8888FF', 'LineWidth', 4, 'DisplayName', '\textbf{Filtered}');
hold off;
title('\boldmath$\mathbf{Cosine \ Signal}$');
xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb. \ units]}$');
legend('show', 'Location', 'best');
xlim([-15, 15]);
ylim([-1.1, 1.1]);
set(gca, 'XColor', 'k');
PlotUtils.styleAxes(gca);
