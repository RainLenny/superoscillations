%% PATH SETUP & CONSTANTS
clear; clc; close all;

%% SETTINGS FOR PARAMETER SWEEP
OmegaTilde_array = linspace(0.001, 0.05, 25); % Array of signal magnitudes to test
num_rng_seeds = 10; % Number of random phases for Rand signal

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

%% BASE SIGNALS SETUP
signal_scaling = 7;

[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, false);

% Flat spectrum signal with same frequencies
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, false);

N_SO = length(angular_freqs_SO);
data_baranov = load('SO_Baranov.mat', 'amps_SO');
amps_baranov = data_baranov.amps_SO * signal_scaling;

%% DATA STORAGE
max_Pe_SO = zeros(length(OmegaTilde_array), 1);
max_Pe_FLAT = zeros(length(OmegaTilde_array), 1);
max_Pe_RAND = zeros(length(OmegaTilde_array), num_rng_seeds);

%% DYNAMICS COMPUTATION
fprintf('Starting parameter sweep over OmegaTilde...\n');

for i = 1:length(OmegaTilde_array)
    fprintf('Computing for OmegaTilde = %.4f (%d/%d)\n', OmegaTilde_array(i), i, length(OmegaTilde_array));
    
    current_params = params;
    current_params.OmegaTilde = OmegaTilde_array(i);
    
    for seed = 1:num_rng_seeds
        [Rand_signal, ~] = generate_rand_phase(angular_freqs_SO, amps_baranov, true, false, seed);
        
        % Normalize signals symbolically by the peak of the first signal (SO)
        [SO_norm, Flat_norm, Rand_norm] = normalize_signals({SO_signal, Flat_signal, Rand_signal}, 'peak');
        
        % For seed == 1, compute SO and Flat
        if seed == 1
            % SO
            [~, rho_SO] = optical_bloch(t_grid, rho_init, nu0, current_params, SO_norm);
            Pe_SO = (1 + rho_SO(:, 3)) / 2;
            max_Pe_SO(i) = max(Pe_SO);
            
            % FLAT
            [~, rho_FLAT] = optical_bloch(t_grid, rho_init, nu0, current_params, Flat_norm);
            Pe_FLAT = (1 + rho_FLAT(:, 3)) / 2;
            max_Pe_FLAT(i) = max(Pe_FLAT);
        end
        
        % RAND (computed for all seeds)
        [~, rho_RAND] = optical_bloch(t_grid, rho_init, nu0, current_params, Rand_norm);
        Pe_RAND = (1 + rho_RAND(:, 3)) / 2;
        max_Pe_RAND(i, seed) = max(Pe_RAND);
    end
end

fprintf('Computation finished.\n');

%% PLOTS
PlotUtils.setupDefaults();

figure('Color', 'w', 'Name', 'Maximal Excitation Comparison');
hold on;

% Plot SO and Flat
plot(OmegaTilde_array, max_Pe_SO, 'LineWidth', 4, 'Color', 'r', 'DisplayName', '\textbf{SO}');
plot(OmegaTilde_array, max_Pe_FLAT, 'LineWidth', 4, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');

% Plot mean and standard deviation for Rand Phase
mean_RAND = mean(max_Pe_RAND, 2);
std_RAND = std(max_Pe_RAND, 0, 2);

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
xlabel('\boldmath$\tilde{\Omega}$ \textbf{(Signal Magnitude)}');
ylabel('\boldmath$\max(P_e)$ \textbf{(Maximal Excitation)}');
legend('show', 'Location', 'northwest');
title('\textbf{Maximal Excitation vs. Signal Magnitude}');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);
