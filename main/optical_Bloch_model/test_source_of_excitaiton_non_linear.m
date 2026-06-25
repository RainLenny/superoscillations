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

%% 1. Define Signals
signal_scaling = 7;

[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, false);
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, false);
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, false);

data_baranov = load('SO_Baranov.mat', 'amps_SO');
amps_baranov = data_baranov.amps_SO * signal_scaling;
[Rand_signal, ~] = generate_rand_phase(angular_freqs_SO, amps_baranov, true, false, 3);

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {}, 'marker', {});

% --- Signal 1: SO ---
sig_configs(1).name = '\textbf{SO}';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;
sig_configs(1).marker = 'o-';

% --- Signal 2: COS ---
sig_configs(2).name = '\boldmath$\mathrm{0.9\omega_0}$';
sig_configs(2).data = Cos_signal;
sig_configs(2).color = 'b';
sig_configs(2).freqs = angular_freqs_COS;
sig_configs(2).marker = 's-';

% --- Signal 3: Flat ---
sig_configs(3).name = '\textbf{Flat}';
sig_configs(3).data = Flat_signal;
sig_configs(3).color = [0, 0.5, 0];
sig_configs(3).freqs = angular_freqs_SO;
sig_configs(3).marker = 'd-';

% --- Signal 4: Rand Phase ---
sig_configs(4).name = '\textbf{Rand Phase}';
sig_configs(4).data = Rand_signal;
sig_configs(4).color = 'm';
sig_configs(4).freqs = angular_freqs_SO;
sig_configs(4).marker = '^-';

% Normalize all signals symbolically by the peak of the first signal
sigs = {sig_configs.data};
[norm_sigs{1:length(sigs)}] = normalize_signals(sigs, 'peak');
for i = 1:length(sig_configs)
    sig_configs(i).data = norm_sigs{i};
end

% Apply defaults (auto-colors)
sig_configs = prepare_signal_config(sig_configs);


%% Dynamics computation
dt = 0.01;
t_grid = (0:dt:T_final)';

% We test different scaling factors lambda
lambdas = logspace(-2, 1, 20); % Range of lambda values

for i = 1:length(sig_configs)
    sig_configs(i).Pe_max = zeros(size(lambdas));
end

for j = 1:length(lambdas)
    lambda = lambdas(j);
    
    for i = 1:length(sig_configs)
        % Scale the signal by lambda
        scaled_signal = @(t) lambda * sig_configs(i).data(t);
        
        [~, rho_out] = optical_bloch(t_grid, rho_init, nu0, params, scaled_signal);
        sig_configs(i).Pe_max(j) = max((1 + rho_out(:, 3)) / 2);
    end

    fprintf('Completed lambda = %.2f\n', lambda);
end


%% PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Amplitude Scaling Log-Log Plot ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test');
hold on;

for i = 1:length(sig_configs)
    loglog(lambdas, sig_configs(i).Pe_max, sig_configs(i).marker, 'LineWidth', 3, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

% Plot reference slopes
ref_lambda = lambdas;
% Align references to the first index of SO (first signal) for visualization
ref_1st_order = sig_configs(1).Pe_max(1) * (ref_lambda / ref_lambda(1)).^2;
ref_3rd_order = sig_configs(1).Pe_max(1) * (ref_lambda / ref_lambda(1)).^6;

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
hold on;

lambda_mid = exp((log(lambdas(1:end-1)) + log(lambdas(2:end))) / 2);

for i = 1:length(sig_configs)
    slope_val = diff(log(sig_configs(i).Pe_max)) ./ diff(log(lambdas));
    semilogx(lambda_mid, slope_val, sig_configs(i).marker, 'LineWidth', 3, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

yline(2, '--k', 'LineWidth', 2, 'DisplayName', '1st order (slope=2)');
yline(6, '-.k', 'LineWidth', 2, 'DisplayName', '3rd order (slope=6)');

grid on;
xlabel('\boldmath$\lambda$ (Drive Amplitude)');
ylabel('\boldmath$d(\log P_e) / d(\log \lambda)$ (Slope)');
title('\textbf{Log-Log Slope vs Amplitude}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);