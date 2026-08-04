%% PATH SETUP & CONSTANTS
clear; clc; close all;

%% SETTINGS FOR PARAMETER SWEEP
OmegaTilde_array = linspace(0.015, 0.025, 50); % Array of signal magnitudes to test
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
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, true);

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {}, 'amps', {}, 'is_rand_phase', {});

% --- Signal 1: SO ---
sig_configs(1).name = '\textbf{SO}';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;
sig_configs(1).amps = amps_SO;
sig_configs(1).is_rand_phase = false;

% % --- Signal 2: Flat ---
% sig_configs(2).name = '\textbf{Flat}';
% sig_configs(2).data = Flat_signal;
% sig_configs(2).color = [0, 0.5, 0];
% sig_configs(2).freqs = angular_freqs_SO;
% sig_configs(2).amps = amps_SO;
% sig_configs(2).is_rand_phase = false;

% --- Signal 3: Rand Phase ---
sig_configs(2).name = '\textbf{Rand Phase}';
sig_configs(2).data = @(t) 0; % Dummy, generated inside loop
sig_configs(2).color = '#008000';
sig_configs(2).freqs = angular_freqs_SO;
sig_configs(2).amps = amps_SO;
sig_configs(2).is_rand_phase = true;

% Apply defaults (auto-colors)
sig_configs = prepare_signal_config(sig_configs);

%% DATA STORAGE
for i = 1:length(sig_configs)
    if sig_configs(i).is_rand_phase
        sig_configs(i).max_Pe = zeros(length(OmegaTilde_array), num_rng_seeds);
        sig_configs(i).Pe_500 = zeros(length(OmegaTilde_array), num_rng_seeds);
    else
        sig_configs(i).max_Pe = zeros(length(OmegaTilde_array), 1);
        sig_configs(i).Pe_500 = zeros(length(OmegaTilde_array), 1);
    end
end

%% DYNAMICS COMPUTATION
fprintf('Starting parameter sweep over OmegaTilde...\n');

for j = 1:length(OmegaTilde_array)
    fprintf('Computing for OmegaTilde = %.4f (%d/%d)\n', OmegaTilde_array(j), j, length(OmegaTilde_array));
    
    
    for seed = 1:num_rng_seeds
        % Generate signals for this seed
        current_data = cell(1, length(sig_configs));
        for i = 1:length(sig_configs)
            if sig_configs(i).is_rand_phase
                [Rand_signal, ~] = generate_rand_phase(sig_configs(i).freqs, sig_configs(i).amps, true, true, seed);
                current_data{i} = Rand_signal;
            else
                current_data{i} = sig_configs(i).data;
            end
        end
        
        % Normalize signals
        [norm_data{1:length(current_data)}] = normalize_signals(current_data, 'energy');
        
        % Simulate
        for i = 1:length(sig_configs)
            if sig_configs(i).is_rand_phase
                [~, Pe] = JC_drive_only(nu0, OmegaTilde_array(j), norm_data{i}, t_grid);
                sig_configs(i).max_Pe(j, seed) = max(Pe);
                sig_configs(i).Pe_500(j, seed) = Pe(idx_t500);
            elseif seed == 1 % only simulate once for non-rand
                [~, Pe] = JC_drive_only(nu0, OmegaTilde_array(j), norm_data{i}, t_grid);
                sig_configs(i).max_Pe(j) = max(Pe);
                sig_configs(i).Pe_500(j) = Pe(idx_t500);
            end
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
                plot(OmegaTilde_array, sig_configs(i).max_Pe(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'DisplayName', '\textbf{Rand Phase (individual)}');
            else
                plot(OmegaTilde_array, sig_configs(i).max_Pe(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end
        % Plot mean
        mean_val = mean(sig_configs(i).max_Pe, 2);
        plot(OmegaTilde_array, mean_val, '--', 'Color', 'k', 'DisplayName', '\textbf{Rand Phase (mean)}');
    else
        plot(OmegaTilde_array, sig_configs(i).max_Pe, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
    end
end

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel(' \textbf{Maximal excitation probability}');
legend('show', 'Location', 'northwest', 'Interpreter', 'latex');
xlim([OmegaTilde_array(1),OmegaTilde_array(end)])



% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);
ax = gca;
ax.XAxis.Exponent = -2;
xtickformat(ax, '$\\mathbf{%g}$');
xticks(linspace(OmegaTilde_array(1), OmegaTilde_array(end), 5));
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
                plot(OmegaTilde_array, sig_configs(i).Pe_500(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'DisplayName', '\textbf{Rand Phase (individual)}');
            else
                plot(OmegaTilde_array, sig_configs(i).Pe_500(:, seed), '-', 'Color', [sig_configs(i).color 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end
        % Plot mean
        mean_val = mean(sig_configs(i).Pe_500, 2);
        plot(OmegaTilde_array, mean_val, '--', 'Color', 'k', 'DisplayName', '\textbf{Rand Phase (Mean)}');
    else
        plot(OmegaTilde_array, sig_configs(i).Pe_500, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
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
xticks(linspace(OmegaTilde_array(1), OmegaTilde_array(end), 5));
xtickangle(0.1);