clear; clc;

%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(mfilename('fullpath'))));
PlotUtils.setupDefaults();

% Folder paths
data_folder = fullfile(fileparts(mfilename('fullpath')), 'data');

%% Load notch filter response
sa_notch_file = fullfile(data_folder, 'sa_notch.csv');
% Read matrix directly (skipping the 34 lines of metadata)
sa_notch_data = readmatrix(sa_notch_file, 'NumHeaderLines', 34);
% Remove last 6 rows as they contain footer info
sa_notch_data(end-5:end, :) = [];

sa_freq = sa_notch_data(:, 1);
sa_dB = sa_notch_data(:, 2);
sa_V = 10.^(sa_dB / 20);

%% Setup variables for looping over notch data
files = {'notch_400khz1.csv', 'notch_500khz1.csv', 'notch_600khz0.csv', 'notch_800khz0.csv'};
freq_labels = {'400 kHz', '500 kHz', '600 kHz', '800 kHz'};
time_limits = [4e-5, 4e-5, 4e-5, 3e-5];

tInc = 1e-9; % Time increment as used in python notebook

for i = 1:length(files)
    %% Load superoscillation data
    filepath = fullfile(data_folder, files{i});
    
    % Use readmatrix, skip 1 header line
    data = readmatrix(filepath, 'NumHeaderLines', 1);
    CH1V = data(:, 1);
    CH2V = data(:, 2);
    
    N = length(CH1V);
    t = (0:N-1)' * tInc;
    
    %% Compute FFT
    fft1 = fft(CH1V);
    fft2 = fft(CH2V);
    freqs = (0:N-1)' / (N * tInc);
    
    fft1_mag = abs(fft1);
    fft2_mag = abs(fft2);
    
    % Normalize
    fft1_norm = 0.9 * fft1_mag / max(fft1_mag);
    fft2_norm = 0.9 * fft2_mag / max(fft2_mag);
    
    %% Plotting
    figure;
    
    % Subplot 1: Time domain
    subplot(2, 1, 1);
    hold on;
    % Python equivalent: notch['CH1V']*300-1 and notch['CH2V']+2.5
    plot(t, CH1V * 300 - 1, '-', 'Color', 'r', 'DisplayName', '\textbf{SO}');
    plot(t, CH2V + 2.5, '-', 'Color', 'b', 'DisplayName', '\textbf{COS}');
    
    title(['\textbf{superoscillation at ', freq_labels{i}, '}']);
    xlabel('\textbf{Time [s]}');
    ylabel('\textbf{Amplitude [V]}');
    xlim([0, time_limits(i)]);
    ylim([-5, 5]);
    legend('show');
    PlotUtils.styleAxes(gca);
    hold off;
    
    % Subplot 2: Frequency domain
    subplot(2, 1, 2);
    hold on;
    % Python equivalent: sa_notch['V']/2.3
    plot(sa_freq, sa_V / 2.3, '-', 'Color', 'k', 'DisplayName', '\textbf{Filter}');
    plot(freqs, fft1_norm, '-o', 'Color', 'r', 'DisplayName', '\textbf{SO}');
    plot(freqs, fft2_norm, '-o', 'Color', 'b', 'DisplayName', '\textbf{COS}');
    
    title(['\textbf{', freq_labels{i}, ' vs. notch frequency response}']);
    xlabel('\textbf{Frequency [Hz]}');
    ylabel('\textbf{Amplitude [V, scaled]}');
    xlim([0, 1e6]);
    ylim([0, 1]);
    legend('show');
    PlotUtils.styleAxes(gca);
    hold off;
end
