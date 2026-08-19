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

legend('show');

% Apply axes styling as in the reference script
PlotUtils.styleAxes(gca);
hold off;
