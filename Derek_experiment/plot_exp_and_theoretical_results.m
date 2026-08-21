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
theory_dB = 20*log10(abs(rlc_data.H_opt(theory_freq))) + rlc_data.P0_final;
theory_V = 10.^(theory_dB / 20);

%% Setup variables for looping over notch data
files = {'notch_400khz1.csv', 'notch_500khz1.csv', 'notch_600khz0.csv', 'notch_800khz0.csv'};
freq_labels = {'400 kHz', '500 kHz', '600 kHz', '800 kHz'};
target_freqs = [400e3, 500e3, 600e3, 800e3]; % Corresponding base frequencies

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
    [~, angular_freqs_SO, amps_SO] = generate_SO_from_dat_file('2p5_cos_Derek', freq_scaling, amp_scaling, false);
    
    % Apply theoretical filter to the SO signal
    f_SO = angular_freqs_SO / (2 * pi);
    H_k = rlc_data.H_opt(f_SO);
    H_k = reshape(H_k, size(f_SO));
    % shift for baseline gain:
    H_k = H_k * 10^(rlc_data.P0_final / 20);
    amps_filtered = amps_SO .* H_k;
    
    [SO_signal_filtered, ~, ~] = generate_signal_base(angular_freqs_SO, amps_filtered, false, false, false);
    SO_theory = SO_signal_filtered(t);
    
    %% Sync and scale SO theory and experimental plots 
    exp_signal = CH1V * 300 - 1; % Scaling from Derek's jupiter notebook
    
    % Align theoretical SO signal
    [r, lags] = xcorr(exp_signal, SO_theory);
    [~, max_idx] = max(r);
    delay_idx = lags(max_idx);
    SO_theory_aligned = circshift(SO_theory, delay_idx);
    

    
    % Scale theoretical time-domain signal to align vertically
    p_so = polyfit(SO_theory_aligned, exp_signal, 1);
    % Add scaling and DC shift:
    SO_theory_time = p_so(1) * SO_theory_aligned + p_so(2);
    
    % Compute theoretical FFT and normalize to fit the experimental measured FFT
    fft_SO_theory = fft(SO_theory);
    fft_SO_theory_mag = abs(fft_SO_theory);
    fft_SO_theory_norm = 0.9 * fft_SO_theory_mag / max(fft_SO_theory_mag);
    
    %% Calculate Zoom Window for Time Domain plot
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
    
    %% Align COS signal to the first peak of the SO in the window
    idx_start = max(1, center_idx - margin);
    idx_end = min(N, center_idx + margin);
    
    % Find peaks in the theoretical SO signal within the window to avoid noise issues
    [SO_pks, SO_locs] = findpeaks(SO_theory_aligned(idx_start:idx_end));
    
    % Filter out small peaks (we want the prominent ones in the superoscillation cluster)
    threshold = 0.5 * max(SO_pks);
    main_SO_locs = SO_locs(SO_pks > threshold);
    
    if ~isempty(main_SO_locs)
        % Get the absolute index of the first prominent peak
        first_SO_peak_idx = idx_start - 1 + main_SO_locs(1);
        
        % Find all true peaks in the blue signal (cosine), ignoring noise ripples
        prominence_thresh = 0.5 * (max(CH2V) - min(CH2V));
        [~, cos_locs] = findpeaks(smoothdata(CH2V, 'gaussian', 15), 'MinPeakProminence', prominence_thresh);
        
        if ~isempty(cos_locs)
            % Find the closest cosine peak to the first SO peak
            [~, closest_idx] = min(abs(cos_locs - first_SO_peak_idx));
            
            % Shift the cosine signal so its peak aligns with the first red peak
            delay_idx_cos = first_SO_peak_idx - cos_locs(closest_idx);
            CH2V_aligned = circshift(CH2V, delay_idx_cos);
        else
            CH2V_aligned = CH2V; % Fallback
        end
    else
        CH2V_aligned = CH2V; % Fallback
    end
    
    %% Plotting
    figure;
    tlo = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Subplot 1: Time domain
    ax1 = nexttile;
    hold on;
    % Plot theoretical SO as dashed so it overlays nicely over experimental (removed from legend)
    h_so_time = plot(t, exp_signal, '-', 'Color', 'r', 'DisplayName', '\textbf{SO}',LineWidth=3);
    h_cos_time = plot(t, CH2V_aligned + 2.5, '-', 'Color', 'b', 'DisplayName', '\textbf{COS}',LineWidth=3);
    plot(t, SO_theory_time, '-', 'Color', [0.5 0 0], 'HandleVisibility', 'off',LineWidth=2);


    % title(['\textbf{superoscillation at ', freq_labels{i}, '}']);
    xlabel('\textbf{Time [s]}');
    ylabel('\textbf{Amplitude [V]}');
    xlim([t_lim_min, t_lim_max]);
    ylim([-2.5, 3.5]);
    PlotUtils.styleAxes(gca);
    hold off;
    
    % Subplot 2: Frequency domain
    ax2 = nexttile;
    hold on;

    % Experimental notch filter response
    h_filter = plot(sa_freq, sa_V / 2.3, '-', 'Color', 'k', 'DisplayName', '\textbf{Filter}');
    % Theoretical notch filter response (dashed, removed from legend)
    plot(theory_freq, theory_V / 2.3, '-', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off',LineWidth=1);
    
    % Experimental SO FFT
    plot(freqs, fft1_norm, '-o', 'Color', 'r', 'DisplayName', '\textbf{SO}',LineWidth=3);
    
    % Experimental COS FFT
    plot(freqs, fft2_norm, '-o', 'Color', 'b', 'DisplayName', '\textbf{COS}',LineWidth=3);
    
    % title(['\textbf{', freq_labels{i}, ' vs. notch frequency response}']);
    xlabel('\textbf{Frequency [Hz]}');
    ylabel('\textbf{Amplitude [V, scaled]}');
    xlim([0, 1e6]);
    ylim([0, 1]);
    
    % Create a single combined legend assigned to the tiledlayout
    lgd_all = legend(ax1, [h_filter, h_so_time, h_cos_time], '\textbf{Filter}', '\textbf{SO}', '\textbf{COS}', 'Orientation', 'horizontal', 'Box', 'off');
    lgd_all.Layout.Tile = 'south';
    
    PlotUtils.styleAxes(gca);
    hold off;
end
