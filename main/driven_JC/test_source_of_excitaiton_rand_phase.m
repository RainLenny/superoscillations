%% PATH SETUP & CONSTANTS
clear; clc; close all;

%% SETTINGS FOR PARAMETER SWEEP
Coupling_array = linspace(0.015, 0.025, 50); % Array of signal magnitudes to test
num_rng_seeds = 40; % Number of random phases for Rand signal

%% CONSTANTS
% Physical and simulation parameters
nu0      = 1.00;

T_final = 500;
dt = 0.01;
t_grid = (0:dt:T_final)';

% Find the index corresponding to t = 500
[~, idx_t500] = min(abs(t_grid - 500));

%% 1. Define Signals
signal_scaling = 1;

[SO_signal, angular_freqs_SO, amps_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, true);

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {}, 'amps', {}, 'is_rand_phase', {});

% --- Signal 1: SO ---
sig_configs(1).name = '\textbf{SO}';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;
sig_configs(1).amps = amps_SO;
sig_configs(1).is_rand_phase = false;

% --- Signal 2: COS ---
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, true);
sig_configs(3).name = '\boldmath$\mathbf{0.9\omega_0}$';
sig_configs(3).data = Cos_signal;
sig_configs(3).color = 'b';
sig_configs(3).freqs = angular_freqs_COS;

% --- Signal 3: Rand Phase ---
sig_configs(2).name = '\textbf{Rand Phase}';
sig_configs(2).data = @(t) 0; % Dummy, generated inside loop
sig_configs(2).color = '#008000';
sig_configs(2).freqs = angular_freqs_SO;
sig_configs(2).amps = amps_SO;
sig_configs(2).is_rand_phase = true;

%% DATA STORAGE
for i = 1:length(sig_configs)
    if sig_configs(i).is_rand_phase
        sig_configs(i).max_Pe = zeros(length(Coupling_array), num_rng_seeds);
        sig_configs(i).Pe_500 = zeros(length(Coupling_array), num_rng_seeds);
    else
        sig_configs(i).max_Pe = zeros(length(Coupling_array), 1);
        sig_configs(i).Pe_500 = zeros(length(Coupling_array), 1);
    end
end

%% PRECOMPUTE AND NORMALIZE SIGNALS
fprintf('Precomputing and normalizing signals...\n');

% Collect all raw signals in a flat list to normalize them together
all_raw_signals = {};
signal_indices = []; % To keep track of which config each signal belongs to

for i = 1:length(sig_configs)
    if sig_configs(i).is_rand_phase
        for seed = 1:num_rng_seeds
            [Rand_signal, ~] = generate_rand_phase(sig_configs(i).freqs, sig_configs(i).amps, true, true, seed);
            all_raw_signals{end+1} = Rand_signal;
            signal_indices(end+1) = i;
        end
    else
        all_raw_signals{end+1} = sig_configs(i).data;
        signal_indices(end+1) = i;
    end
end

% Normalize all signals relative to the first signal (SO)
norm_results = cell(1, length(all_raw_signals));
[norm_results{:}] = normalize_signals(all_raw_signals, 'energy');

% Distribute the normalized signals back to sig_configs
for i = 1:length(sig_configs)
    sig_configs(i).norm_signals = norm_results(signal_indices == i);
end

%% DYNAMICS COMPUTATION
fprintf('Starting parameter sweep over OmegaTilde...\n');

for j = 1:length(Coupling_array)
    fprintf('Computing for OmegaTilde = %.4f (%d/%d)\n', Coupling_array(j), j, length(Coupling_array));
    
    for i = 1:length(sig_configs)
        if sig_configs(i).is_rand_phase
            for seed = 1:num_rng_seeds
                [~, Pe] = JC_drive_only(nu0, Coupling_array(j), sig_configs(i).norm_signals{seed}, t_grid);
                sig_configs(i).max_Pe(j, seed) = max(Pe);
                sig_configs(i).Pe_500(j, seed) = Pe(idx_t500);
            end
        else
            % Simulate only once for non-rand
            [~, Pe] = JC_drive_only(nu0, Coupling_array(j), sig_configs(i).norm_signals{1}, t_grid);
            sig_configs(i).max_Pe(j) = max(Pe);
            sig_configs(i).Pe_500(j) = Pe(idx_t500);
        end
    end
end

fprintf('Computation finished.\n');

%% PLOTS
PlotUtils.setupDefaults();

% =========================================================
% PLOT 1: Maximal Excitation
% =========================================================
figure('Color', 'w', 'Name', 'Maximal Excitation Comparison');
hold on;

for i = 1:length(sig_configs)
    if sig_configs(i).is_rand_phase
        % Plot individual random seeds as scattered lines with transparency
        for seed = 1:num_rng_seeds
            if seed == 1
                plot(Coupling_array, sig_configs(i).max_Pe(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'DisplayName', '\textbf{Rand Phase (individual)}');
            else
                plot(Coupling_array, sig_configs(i).max_Pe(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end
        % Plot mean
        mean_val = mean(sig_configs(i).max_Pe, 2);
        plot(Coupling_array, mean_val, '--', 'Color', 'k', 'DisplayName', '\textbf{Rand Phase (mean)}');
    else
        plot(Coupling_array, sig_configs(i).max_Pe, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
    end
end

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel(' \textbf{Maximal excitation probability}');
legend('show', 'Location', 'northwest', 'Interpreter', 'latex');
xlim([Coupling_array(1),Coupling_array(end)])



% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);
ax = gca;
ax.XAxis.Exponent = -2;
xtickformat(ax, '$\\mathbf{%g}$');
xticks(linspace(Coupling_array(1), Coupling_array(end), 5));
xtickangle(0.1);
ylim([0.2,1])

% =========================================================
% PLOT 2: Excitation at t = 500
% =========================================================
figure('Color', 'w', 'Name', 'Excitation at t=500 Comparison');
hold on;

for i = 1:length(sig_configs)
    if sig_configs(i).is_rand_phase
        % Plot individual random seeds as scattered lines with transparency
        for seed = 1:num_rng_seeds
            if seed == 1
                plot(Coupling_array, sig_configs(i).Pe_500(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'DisplayName', '\textbf{Rand Phase (individual)}');
            else
                plot(Coupling_array, sig_configs(i).Pe_500(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end
        % Plot mean
        mean_val = mean(sig_configs(i).Pe_500, 2);
        plot(Coupling_array, mean_val, '--', 'Color', 'k', 'DisplayName', '\textbf{Rand Phase (Mean)}');
    else
        plot(Coupling_array, sig_configs(i).Pe_500, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
    end
end

grid on;
xlabel('\textbf{Signal Magnitude}');
ylabel('\boldmath$P_e(t=500)$ \textbf{(Excitation at t=500)}');
legend('show', 'Location', 'northwest', 'Interpreter', 'latex');
title('\textbf{Excitation at $t=500$ vs. Signal Magnitude}');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);
ax = gca;
ax.XAxis.Exponent = -2;
xtickformat(ax, '$\\mathbf{%g}$');
xticks(linspace(Coupling_array(1), Coupling_array(end), 5));
xtickangle(0.1);