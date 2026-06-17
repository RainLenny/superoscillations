%% PATH SETUP & CONSTANTS
clear; clc; close all;
% Add all project subfolders to search path
addpath(genpath(fileparts(mfilename('fullpath'))));

%% CONSTANTS

% Time and envelope parameters
t_0 = 250; 
T   = 100;

% Physical and simulation parameters
nu0      = 1.00;
rho_init = [0; 0; -1]; %Inital conditions of the TLS

params.T1         = 500;
params.T2         = 300;
params.rho30      = -1.0; %Value to which the system relaxes to
params.OmegaTilde = 0.01; %Field-TLS coupling 

T_final = 600;

%% SIGNALS:
signal_scaling = 7;

[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, false);

[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, false);

% Flat spectrum signal with same frequencies
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, false);

% Random phase signal (keeping Baranov amplitudes but with random phase, real-valued)
data_baranov = load('SO_Baranov.mat', 'amps_SO');
amps_baranov = data_baranov.amps_SO * signal_scaling;
[Rand_signal, ~] = generate_rand_phase(angular_freqs_SO, amps_baranov, true, false, 1);

% Normalize all signals symbolically by the peak of the first signal
[SO_signal, Cos_signal, Flat_signal, Rand_signal] = normalize_signals({SO_signal, Cos_signal, Flat_signal, Rand_signal}, 'peak');

%% Dynamics computation
dt = 0.01;
t_grid = (0:dt:T_final)';

% We test different scaling factors lambda
lambdas = logspace(-2, 1, 20); % Range of lambda values
Pe_max_SO   = zeros(size(lambdas));
Pe_max_COS  = zeros(size(lambdas));
Pe_max_FLAT = zeros(size(lambdas));
Pe_max_RAND = zeros(size(lambdas));

for i = 1:length(lambdas)
    lambda = lambdas(i);
    
    % Scale the signals by lambda
    scaled_SO   = @(t) lambda * SO_signal(t);
    scaled_COS  = @(t) lambda * Cos_signal(t);
    scaled_FLAT = @(t) lambda * Flat_signal(t);
    scaled_RAND = @(t) lambda * Rand_signal(t);
    
    [~, rho_SO] = optical_bloch(t_grid, rho_init, nu0, params, scaled_SO);
    Pe_max_SO(i) = max((1 + rho_SO(:, 3)) / 2);
    
    [~, rho_COS] = optical_bloch(t_grid, rho_init, nu0, params, scaled_COS);
    Pe_max_COS(i) = max((1 + rho_COS(:, 3)) / 2);
    
    [~, rho_FLAT] = optical_bloch(t_grid, rho_init, nu0, params, scaled_FLAT);
    Pe_max_FLAT(i) = max((1 + rho_FLAT(:, 3)) / 2);
    
    [~, rho_RAND] = optical_bloch(t_grid, rho_init, nu0, params, scaled_RAND);
    Pe_max_RAND(i) = max((1 + rho_RAND(:, 3)) / 2);

    fprintf('Completed lambda = %.2f\n', lambda);

end


%% PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Amplitude Scaling Log-Log Plot ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test');
loglog(lambdas, Pe_max_SO, 'o-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', 'red', 'DisplayName', '\textbf{SO}');
hold on;
loglog(lambdas, Pe_max_COS, 's-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', 'blue', 'DisplayName', '\boldmath$\mathrm{0.9\omega_0}$');
loglog(lambdas, Pe_max_FLAT, 'd-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');
loglog(lambdas, Pe_max_RAND, '^-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', 'm', 'DisplayName', '\textbf{Rand Phase}');

% Plot reference slopes
ref_lambda = lambdas;
% Align references to the first index of SO for visualization
ref_1st_order = Pe_max_SO(1) * (ref_lambda / ref_lambda(1)).^2;
ref_3rd_order = Pe_max_SO(1) * (ref_lambda / ref_lambda(1)).^6;

loglog(ref_lambda, ref_1st_order, '--k', 'LineWidth', 2, 'DisplayName', '\boldmath$\propto \lambda^2$ (1st order)');
loglog(ref_lambda, ref_3rd_order, '-.k', 'LineWidth', 2, 'DisplayName', '\boldmath$\propto \lambda^6$ (3rd order)');

grid on;
xlabel('\boldmath$\lambda$ (Drive Amplitude)');
ylabel('\boldmath$\max(P_e)$ (Max Excitation Probability)');

ylim([-inf, 1])

title('\textbf{Amplitude Scaling Test: Log-Log Plot}');
legend('show', 'Location', 'best');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 2: Log-Log Slopes ---
figure('Color', 'w', 'Name', 'Log-Log Slopes');
slope_SO = diff(log(Pe_max_SO)) ./ diff(log(lambdas));
slope_COS = diff(log(Pe_max_COS)) ./ diff(log(lambdas));
slope_FLAT = diff(log(Pe_max_FLAT)) ./ diff(log(lambdas));
slope_RAND = diff(log(Pe_max_RAND)) ./ diff(log(lambdas));

lambda_mid = exp((log(lambdas(1:end-1)) + log(lambdas(2:end))) / 2);

semilogx(lambda_mid, slope_SO, 'o-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', 'red', 'DisplayName', '\textbf{SO}');
hold on;
semilogx(lambda_mid, slope_COS, 's-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', 'blue', 'DisplayName', '\boldmath$\mathrm{0.9\omega_0}$');
semilogx(lambda_mid, slope_FLAT, 'd-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');
semilogx(lambda_mid, slope_RAND, '^-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', 'm', 'DisplayName', '\textbf{Rand Phase}');

yline(2, '--k', 'LineWidth', 2, 'DisplayName', '1st order (slope=2)');
yline(6, '-.k', 'LineWidth', 2, 'DisplayName', '3rd order (slope=6)');

grid on;
xlabel('\boldmath$\lambda$ (Drive Amplitude)');
ylabel('\boldmath$d(\log P_e) / d(\log \lambda)$ (Slope)');
title('\textbf{Log-Log Slope vs Amplitude}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);