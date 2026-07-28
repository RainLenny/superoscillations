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

[SO_signal, angular_freqs_SO, amps_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, false);
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, false);
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, false);

[Rand_signal, ~] = generate_rand_phase(angular_freqs_SO, amps_SO, true, false, 3);

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {}, 'marker', {});

% --- Signal 1: SO ---
sig_configs(1).name = '\textbf{SO}';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;
sig_configs(1).marker = '-';

% --- Signal 2: COS ---
sig_configs(2).name = '\boldmath$\mathbf{0.9\omega_0}$';
sig_configs(2).data = Cos_signal;
sig_configs(2).color = 'b';
sig_configs(2).freqs = angular_freqs_COS;
sig_configs(2).marker = '-';

% % --- Signal 3: Flat ---
% sig_configs(3).name = '\textbf{Flat}';
% sig_configs(3).data = Flat_signal;
% sig_configs(3).color = [0, 0.5, 0];
% sig_configs(3).freqs = angular_freqs_SO;
% sig_configs(3).marker = 'd-';
% 
% % --- Signal 4: Rand Phase ---
% sig_configs(4).name = '\textbf{Rand Phase}';
% sig_configs(4).data = Rand_signal;
% sig_configs(4).color = 'm';
% sig_configs(4).freqs = angular_freqs_SO;
% sig_configs(4).marker = '^-';

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
lambdas = logspace(-1, 0, 100); % Range of lambda values

for i = 1:length(sig_configs)
    sig_configs(i).Pe_max = zeros(size(lambdas));
    sig_configs(i).rho1_max = zeros(size(lambdas));
    sig_configs(i).rho2_max = zeros(size(lambdas));
end

for j = 1:length(lambdas)
    lambda = lambdas(j);
    
    for i = 1:length(sig_configs)
        % Scale the signal by lambda
        scaled_signal = @(t) lambda * sig_configs(i).data(t);
        
        [~, rho_out] = optical_bloch(t_grid, rho_init, nu0, params, scaled_signal);
        sig_configs(i).Pe_max(j) = max((1 + rho_out(:, 3)) / 2);
        sig_configs(i).rho1_max(j) = max(abs(rho_out(:, 1)));
        sig_configs(i).rho2_max(j) = max(abs(rho_out(:, 2)));
    end

    fprintf('Completed lambda = %.2f\n', lambda);
end


%% PLOTS
PlotUtils.setupDefaults();

% Set default interpreters to LaTeX to ensure custom formatting functions flawlessly
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% --- Figure 1: Amplitude Scaling Log-Log Plot ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test');
hold on;

for i = 1:length(sig_configs)
    loglog(lambdas, sig_configs(i).Pe_max, sig_configs(i).marker, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

% Plot reference slopes
ref_lambda = lambdas;
% Align references to the first index of SO (first signal) for visualization
ref_1st_order = sig_configs(1).Pe_max(1) * (ref_lambda / ref_lambda(1)).^2;
ref_3rd_order = sig_configs(1).Pe_max(1) * (ref_lambda / ref_lambda(1)).^6;

loglog(ref_lambda, ref_1st_order, '--k', 'DisplayName', '\boldmath$\propto \lambda^2$ \textbf{(1st order)}');
loglog(ref_lambda, ref_3rd_order, '-.k', 'DisplayName', '\boldmath$\propto \lambda^6$ \textbf{(3rd order)}');

grid on;
xlabel('\boldmath$\lambda$ \textbf{(Drive Amplitude)}');
ylabel('\boldmath$\max(P_e)$ \textbf{(Max Excitation Probability)}');

ylim([-inf, 1])

legend('show', 'Location', 'best');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold'); % Enforce bold ticks/scales post-styling


% --- Figure 2: Log-Log Slopes ---
figure('Color', 'w', 'Name', 'Log-Log Slopes');
hold on;

lambda_mid = exp((log(lambdas(1:end-1)) + log(lambdas(2:end))) / 2);

for i = 1:length(sig_configs)
    slope_val = diff(log(sig_configs(i).Pe_max)) ./ diff(log(lambdas));
    semilogx(lambda_mid, slope_val, sig_configs(i).marker, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

yline(2, '--k', 'DisplayName', '\textbf{1st order (slope=2)}');
yline(6, '-.k', 'DisplayName', '\textbf{3rd order (slope=6)}');

grid on;
xlabel('\boldmath$\lambda$ \textbf{(Drive Amplitude)}');
ylabel('\boldmath$d(\log P_e) / d(\log \lambda)$ \textbf{(Slope)}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'FontWeight', 'bold');


% --- Figure 3: Amplitude Scaling Log-Log Plot (rho_1) ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test (rho_1)');
hold on;

for i = 1:length(sig_configs)
    loglog(lambdas, sig_configs(i).rho1_max, sig_configs(i).marker, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

% Plot reference slopes
ref_lambda = lambdas;
ref_1st_order_rho = sig_configs(1).rho1_max(1) * (ref_lambda / ref_lambda(1)).^1;
ref_3rd_order_rho = sig_configs(1).rho1_max(1) * (ref_lambda / ref_lambda(1)).^3;

loglog(ref_lambda, ref_1st_order_rho, '--k', 'DisplayName', '\boldmath$\propto \lambda^1$ \textbf{(1st order)}');
loglog(ref_lambda, ref_3rd_order_rho, '-.k', 'DisplayName', '\boldmath$\propto \lambda^3$ \textbf{(3rd order)}');

grid on;
xlabel('\boldmath$\lambda$ \textbf{(Drive Amplitude)}');
ylabel('\boldmath$\max(|\rho_1|)$');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold');


% --- Figure 4: Log-Log Slopes (rho_1) ---
figure('Color', 'w', 'Name', 'Log-Log Slopes (rho_1)');
hold on;

for i = 1:length(sig_configs)
    slope_val = diff(log(sig_configs(i).rho1_max)) ./ diff(log(lambdas));
    semilogx(lambda_mid, slope_val, sig_configs(i).marker, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

yline(1, '--k', 'DisplayName', '\textbf{1st order (slope=1)}');
yline(3, '-.k', 'DisplayName', ' \textbf{3rd order (slope=3)}');

grid on;
xlabel('\boldmath$\lambda$ \textbf{(Drive Amplitude)}');
ylabel('\boldmath$d(\log |\rho_1|) / d(\log \lambda)$ \textbf{(Slope)}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'FontWeight', 'bold');


% --- Figure 5: Amplitude Scaling Log-Log Plot (rho_2) ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test (rho_2)');
hold on;

for i = 1:length(sig_configs)
    loglog(lambdas, sig_configs(i).rho2_max, sig_configs(i).marker, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

% Plot reference slopes
ref_1st_order_rho2 = sig_configs(1).rho2_max(1) * (ref_lambda / ref_lambda(1)).^1;
ref_3rd_order_rho2 = sig_configs(1).rho2_max(1) * (ref_lambda / ref_lambda(1)).^3;

loglog(ref_lambda, ref_1st_order_rho2, '--k', 'DisplayName', '\boldmath$\propto \lambda^1$ \textbf{(1st order)}');
loglog(ref_lambda, ref_3rd_order_rho2, '-.k', 'DisplayName', '\boldmath$\propto \lambda^3$ \textbf{(3rd order)}');

grid on;
xlabel('\boldmath$\lambda$ \textbf{(Drive Amplitude)}');
ylabel('\boldmath$\max(|\rho_2|)$');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold');


% --- Figure 6: Log-Log Slopes (rho_2) ---
figure('Color', 'w', 'Name', 'Log-Log Slopes (rho_2)');
hold on;

for i = 1:length(sig_configs)
    slope_val = diff(log(sig_configs(i).rho2_max)) ./ diff(log(lambdas));
    semilogx(lambda_mid, slope_val, sig_configs(i).marker, 'MarkerSize', 8, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

yline(1, '--', 'DisplayName', '\textbf{1st order (slope=1)}');
yline(3, '-.', 'DisplayName', '\textbf{3rd order (slope=3)}');

grid on;
xlabel('\boldmath$\lambda$ \textbf{(Drive Amplitude)}');
ylabel('\boldmath$d(\log |\rho_2|) / d(\log \lambda)$ \textbf{(Slope)}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'FontWeight', 'bold');