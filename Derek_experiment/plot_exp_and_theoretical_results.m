clear; clc;

%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(mfilename('fullpath'))));
PlotUtils.setupDefaults();

% Folder paths
data_folder = fullfile(fileparts(mfilename('fullpath')), 'data');

%% Load notch filter response (experimental)
sa_notch_file = fullfile(data_folder, 'sa_notch.csv');
% Read matrix directly (skipping the 34 lines of metadata)
sa_notch_data = readmatrix(sa_notch_file, 'NumHeaderLines', 34);
% Remove last 6 rows as they contain footer info
sa_notch_data(end-5:end, :) = [];

sa_freq = sa_notch_data(:, 1);
sa_dB = sa_notch_data(:, 2);
sa_V = 10.^(sa_dB / 20);

%% Load theoretical notch filter response
rlc_data = load(fullfile(fileparts(mfilename('fullpath')), 'optimized_RLC.mat'));
theory_freq = rlc_data.freq_full;
theory_dB = 20*log10(abs(rlc_data.H_opt)) + rlc_data.P0_final;
theory_V = 10.^(theory_dB / 20);

%% Setup variables for looping over notch data
files = {'notch_400khz1.csv', 'notch_500khz1.csv', 'notch_600khz0.csv', 'notch_800khz0.csv'};
freq_labels = {'400 kHz', '500 kHz', '600 kHz', '800 kHz'};
target_freqs = [400e3, 500e3, 600e3, 800e3]; % Corresponding base frequencies
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
    
    %% Compute experimental FFT
    fft1 = fft(CH1V);
    fft2 = fft(CH2V);
    freqs = (0:N-1)' / (N * tInc);
    
    fft1_mag = abs(fft1);
    fft2_mag = abs(fft2);
    
    % Normalize experimental FFT
    fft1_norm = 0.9 * fft1_mag / max(fft1_mag);
    fft2_norm = 0.9 * fft2_mag / max(fft2_mag);
    
    %% Generate theoretical SO signal
    % target_freqs correspond to local frequency scale (e.g. 1.0 -> target_freq)
    freq_scaling = 2 * pi * target_freqs(i); 
    amp_scaling = 1;
    [SO_signal, ~, ~] = generate_SO_from_dat_file('2p5_cos_Derek', freq_scaling, amp_scaling, false);
    
    SO_theory = SO_signal(t);
    
    % Calculate the time delay seen in the time plot to sync up the plots
    exp_signal = CH1V * 300 - 1;
    [r, lags] = xcorr(exp_signal, SO_theory);
    [~, max_idx] = max(r);
    delay_idx = lags(max_idx);
    
    % Shift the theoretical signal to sync with the experimental one
    SO_theory_aligned = circshift(SO_theory, delay_idx);
    
    % Scale theoretical time-domain signal to align vertically
    p_so = polyfit(SO_theory_aligned, exp_signal, 1);
    SO_theory_time = p_so(1) * SO_theory_aligned + p_so(2);
    
    % Compute theoretical FFT and normalize to fit the experimental measured FFT
    fft_SO_theory = fft(SO_theory);
    fft_SO_theory_mag = abs(fft_SO_theory);
    fft_SO_theory_norm = 0.9 * fft_SO_theory_mag / max(fft_SO_theory_mag);
    
    %% Calculate Zoom Window for Time Domain
    % The signal has an envelope period T_env determined by the frequency step.
    % Base frequencies are [0.3, 0.4, 0.5, 0.6], step is 0.1.
    % Scaled frequency step is 0.1 * target_freqs(i).
    T_env = 1 / (0.1 * target_freqs(i));
    N_env = round(T_env / tInc);
    
    % Find a main lobe to center the plot
    [~, global_max_idx] = max(abs(SO_theory_aligned));
    
    % Main lobes repeat every N_env. Find all main lobes in the data.
    k_min = -ceil(global_max_idx / N_env);
    k_max = ceil((N - global_max_idx) / N_env);
    main_lobes = global_max_idx + (k_min:k_max) * N_env;
    
    % We want a window of width 1.5 * T_env centered on a main lobe,
    % which will capture exactly two superoscillation clusters (one on each side).
    margin = round(0.75 * N_env);
    valid_lobes = main_lobes((main_lobes - margin > 0) & (main_lobes + margin <= N));
    
    if ~isempty(valid_lobes)
        center_idx = valid_lobes(1); % Pick the first fully visible main lobe
    else
        center_idx = round(N/2); % Fallback
    end
    
    t_lim_min = t(max(1, center_idx - margin));
    t_lim_max = t(min(N, center_idx + margin));
    
    %% Plotting
    figure;
    
    % Subplot 1: Time domain
    subplot(2, 1, 1);
    hold on;
    % Plot theoretical SO as dashed so it overlays nicely over experimental (removed from legend)
    plot(t, exp_signal, '-', 'Color', 'r', 'DisplayName', '\textbf{SO}');
    plot(t, CH2V + 2.5, '-', 'Color', 'b', 'DisplayName', '\textbf{COS}');
    plot(t, SO_theory_time, ':', 'Color', [0.5 0 0], 'HandleVisibility', 'off');


    title(['\textbf{superoscillation at ', freq_labels{i}, '}']);
    xlabel('\textbf{Time [s]}');
    ylabel('\textbf{Amplitude [V]}');
    xlim([t_lim_min, t_lim_max]);
    ylim([-5, 5]);
    legend('show');
    PlotUtils.styleAxes(gca);
    hold off;
    
    % Subplot 2: Frequency domain
    subplot(2, 1, 2);
    hold on;

    % Experimental notch filter response
    plot(sa_freq, sa_V / 2.3, '-', 'Color', 'k', 'DisplayName', '\textbf{Filter}');
    % Theoretical notch filter response (dashed, removed from legend)
    plot(theory_freq, theory_V / 2.3, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    
    % Experimental SO FFT
    plot(freqs, fft1_norm, '-o', 'Color', 'r', 'DisplayName', '\textbf{SO}');
    
    % Experimental COS FFT
    plot(freqs, fft2_norm, '-o', 'Color', 'b', 'DisplayName', '\textbf{COS}');
    
    % title(['\textbf{', freq_labels{i}, ' vs. notch frequency response}']);
    xlabel('\textbf{Frequency [Hz]}');
    ylabel('\textbf{Amplitude [V, scaled]}');
    xlim([0, 1e6]);
    ylim([0, 1]);
    legend('show');
    PlotUtils.styleAxes(gca);
    hold off;
end
