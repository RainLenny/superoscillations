%% Setup Parameters
clear; clc; close all;

%% Global constants
%Time Span
t_span = [0 600];

% 2. Initial Conditions
rho_init = [0; 0; -1]; 

% 3. System Constants
params.T1      = 500;  % Longitudinal relaxation
params.T2      = 300;   % Transverse relaxation
params.rho30   = -1.0;   % Equilibrium for rho3
params.OmegaTilde = 0.01;

nu0 = 1.01;

%% Signal
[f_signal, ~, ~] = generate_signals_Baranov_2014(1,7);
%% Run Solver
[t, rho] = solve_bloch_odes(t_span, rho_init, nu0, params, f_signal);

%% Visualize Results
figure('Color', 'w');
% plot(t, rho(:,1), 'r', 'LineWidth', 1.5, 'DisplayName', '\rho_1'); hold on;
% plot(t, rho(:,2), 'g', 'LineWidth', 1.5, 'DisplayName', '\rho_2');
% yline(params.rho30, '--k', 'DisplayName', '\rho_{30} (Eq.)');
plot(t, rho(:,3), 'b', 'LineWidth', 1.5, 'DisplayName', '\rho_3');


xlabel('Time (t)');
ylabel('Amplitude');
title('Dynamics of Coupled ODEs');
legend('Location', 'best');
grid on;