clear; clc;
%% Initialization
% Add project root and all subfolders to search path
% This script is in Derek_experiment/plot/
% Two fileparts up gets us to Derek_experiment, three gets us to root ('for chat')
addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
PlotUtils.setupDefaults();

%% Load Data
% Get the path to sa_notch.csv
dataFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'sa_notch.csv');

% Read data, skipping the 34 header lines
% readmatrix will return NaN for the footer rows which contain strings
data = readmatrix(dataFile, 'NumHeaderLines', 34);

% Filter out any NaN rows (the footer lines)
valid_idx = ~isnan(data(:,1)) & ~isnan(data(:,2));
freq_Hz = data(valid_idx, 1);
amp_dBm = data(valid_idx, 2);

%% Plot the Filter
figure;
hold on;

% Using the same conventions: 'k' color for filter, LaTeX formatted DisplayName
plot(freq_Hz, amp_dBm, '-', 'Color', 'k', 'DisplayName', '\textbf{Filter}');

% LaTeX formatted labels as in the reference script
xlabel('\boldmath$\mathbf{Frequency \ [Hz]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [dBm]}$');
ylim([-5, 10]);

legend('show');

% Apply axes styling as in the reference script
PlotUtils.styleAxes(gca);
hold off;

%% Load Theoretical Data
rlc_data = load(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'optimized_RLC.mat'));
theory_freq = rlc_data.freq_full;
% Compute theoretical amplitude in dB (to match the dBm of the experiment)
theory_dB = 20*log10(abs(rlc_data.H_opt(theory_freq))) + rlc_data.P0_final;
% Compute theoretical phase
theory_phase = angle(rlc_data.H_opt(theory_freq));

%% Plot Theoretical Filter and Phase
figure;
tlo = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Magnitude
nexttile;
hold on;
% Plot theoretical filter
plot(theory_freq, theory_dB, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', '\textbf{Theoretical Filter}');
% Plot experimental filter (sa_notch) on top
plot(freq_Hz, amp_dBm, '-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', '\textbf{Experimental Filter}');

xlabel('\boldmath$\mathbf{Frequency \ [Hz]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [dBm]}$');
ylim([-5, 10]);
legend('show');

PlotUtils.styleAxes(gca);
hold off;

% Subplot 2: Phase
nexttile;
hold on;
plot(theory_freq, theory_phase, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', '\textbf{Theoretical Phase}');

xlabel('\boldmath$\mathbf{Frequency \ [Hz]}$');
ylabel('\boldmath$\mathbf{Phase \ [rad]}$');
legend('show');

PlotUtils.styleAxes(gca);
hold off;
