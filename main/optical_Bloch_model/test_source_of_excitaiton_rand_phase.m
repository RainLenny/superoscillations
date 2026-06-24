%% PATH SETUP & CONSTANTS
clear; clc; close all;

%% SETTINGS FOR PARAMETER SWEEP
OmegaTilde_array = linspace(0.005, 0.02, 50); % Array of signal magnitudes to test
num_rng_seeds = 30; % Number of random phases for Rand signal

%% CONSTANTS
% Physical and simulation parameters
nu0      = 1.00;
rho_init = [0; 0; -1]; % Inital conditions of the TLS
params.T1         = 500;
params.T2         = 300;
params.rho30      = -1.0; % Value to which the system relaxes to

T_final = 600;
dt = 0.01;
t_grid = (0:dt:T_final)';

% Find the index corresponding to t = 500
[~, idx_t500] = min(abs(t_grid - 500));

%% BASE SIGNALS SETUP
signal_scaling = 7;

[SO_signal, angular_freqs_SO, amps_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, false);

% Flat spectrum signal with same frequencies
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, false);


%% DATA STORAGE
% Storage for Maximal Excitation
max_Pe_SO = zeros(length(OmegaTilde_array), 1);
max_Pe_FLAT = zeros(length(OmegaTilde_array), 1);
max_Pe_RAND = zeros(length(OmegaTilde_array), num_rng_seeds);

% Storage for Excitation at t = 500
Pe_500_SO = zeros(length(OmegaTilde_array), 1);
Pe_500_FLAT = zeros(length(OmegaTilde_array), 1);
Pe_500_RAND = zeros(length(OmegaTilde_array), num_rng_seeds);

%% DYNAMICS COMPUTATION
fprintf('Starting parameter sweep over OmegaTilde...\n');

for i = 1:length(OmegaTilde_array)
    fprintf('Computing for OmegaTilde = %.4f (%d/%d)\n', OmegaTilde_array(i), i, length(OmegaTilde_array));
    
    current_params = params;
    current_params.OmegaTilde = OmegaTilde_array(i);
    
    for seed = 1:num_rng_seeds
        [Rand_signal, ~] = generate_rand_phase(angular_freqs_SO, amps_SO, true, false, seed);
        
        % Normalize signals symbolically by the peak of the first signal (SO)
        [SO_norm, Flat_norm, Rand_norm] = normalize_signals({SO_signal, Flat_signal, Rand_signal}, 'peak');
        
        % For seed == 1, compute SO and Flat
        if seed == 1
            % SO
            [~, rho_SO] = optical_bloch(t_grid, rho_init, nu0, current_params, SO_norm);
            Pe_SO = (1 + rho_SO(:, 3)) / 2;
            max_Pe_SO(i) = max(Pe_SO);
            Pe_500_SO(i) = Pe_SO(idx_t500); % Store excitation at t = 500
            
            % FLAT
            [~, rho_FLAT] = optical_bloch(t_grid, rho_init, nu0, current_params, Flat_norm);
            Pe_FLAT = (1 + rho_FLAT(:, 3)) / 2;
            max_Pe_FLAT(i) = max(Pe_FLAT);
            Pe_500_FLAT(i) = Pe_FLAT(idx_t500); % Store excitation at t = 500
        end
        
        % RAND (computed for all seeds)
        [~, rho_RAND] = optical_bloch(t_grid, rho_init, nu0, current_params, Rand_norm);
        Pe_RAND = (1 + rho_RAND(:, 3)) / 2;
        max_Pe_RAND(i, seed) = max(Pe_RAND);
        Pe_500_RAND(i, seed) = Pe_RAND(idx_t500); % Store excitation at t = 500
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

% Plot SO and Flat
plot(OmegaTilde_array, max_Pe_SO, 'LineWidth', 4, 'Color', 'r', 'DisplayName', '\textbf{SO}');
plot(OmegaTilde_array, max_Pe_FLAT, 'LineWidth', 4, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');

% Plot mean and standard deviation for Rand Phase
mean_RAND = mean(max_Pe_RAND, 2);

% Plot individual random seeds as scattered lines with transparency
for seed = 1:num_rng_seeds
    if seed == 1
        plot(OmegaTilde_array, max_Pe_RAND(:, seed), '-', 'Color', [1 0 1 0.2], 'LineWidth', 1.5, 'DisplayName', '\textbf{Rand (individual)}');
    else
        plot(OmegaTilde_array, max_Pe_RAND(:, seed), '-', 'Color', [1 0 1 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

% Plot mean of Rand Phase
plot(OmegaTilde_array, mean_RAND, '--', 'Color', 'm', 'LineWidth', 4, 'DisplayName', '\textbf{Rand Phase (Mean)}');

grid on;
xlabel('\textbf{Signal Magnitude}');
ylabel('\boldmath$\max(P_e)$ \textbf{(Maximal Excitation)}');
legend('show', 'Location', 'northwest');
title('\textbf{Maximal Excitation vs. Signal Magnitude}');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% =========================================================
% PLOT 2: Excitation at t = 500
% =========================================================
figure('Color', 'w', 'Name', 'Excitation at t=500 Comparison');
hold on;

% Plot SO and Flat
plot(OmegaTilde_array, Pe_500_SO, 'LineWidth', 4, 'Color', 'r', 'DisplayName', '\textbf{SO}');
plot(OmegaTilde_array, Pe_500_FLAT, 'LineWidth', 4, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');

% Plot mean for Rand Phase at t=500
mean_RAND_500 = mean(Pe_500_RAND, 2);

% Plot individual random seeds as scattered lines with transparency
for seed = 1:num_rng_seeds
    if seed == 1
        plot(OmegaTilde_array, Pe_500_RAND(:, seed), '-', 'Color', [1 0 1 0.2], 'LineWidth', 1.5, 'DisplayName', '\textbf{Rand (individual)}');
    else
        plot(OmegaTilde_array, Pe_500_RAND(:, seed), '-', 'Color', [1 0 1 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

% Plot mean of Rand Phase
plot(OmegaTilde_array, mean_RAND_500, '--', 'Color', 'm', 'LineWidth', 4, 'DisplayName', '\textbf{Rand Phase (Mean)}');

grid on;
xlabel('\textbf{Signal Magnitude}');
ylabel('\boldmath$P_e(t=500)$ \textbf{(Excitation at t=500)}');
legend('show', 'Location', 'northwest');
title('\textbf{Excitation at $t=500$ vs. Signal Magnitude}');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);