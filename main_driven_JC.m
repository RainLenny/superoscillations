%% PATH SETUP & CONSTANTS
clear; clc;
% Add all project subfolders to search path
addpath(genpath(fileparts(mfilename('fullpath'))));

%% CONSTANTS

% Time and envelope parameters
t_0 = 250; 
T   = 100;

% Physical and simulation parameters
nu0    = 1;
J_drive = 6;
J_fluc  = 0.003;

T_final = 500;
nmax    = 0;

% Control flags
do_err_est = 0;

%% SIGNALS:
signal_scaling = 1.3;

[SO_signal, angular_freqs_SO] = generate_SO_Baranov(1,signal_scaling,true,true);

[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9,signal_scaling,true,true);

% Normalize both signals symbolically by the peak of the first signal
[SO_signal, Cos_signal] = normalize_signals({SO_signal, Cos_signal}, 'peak');

%% Dyamics computation
% SO
[tgrid_SO, Pe_SO, nk_SO, ntot_SO, eps_trunc_SO] = multimode_JC_driven_mode_photon_cap(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, SO_signal, T_final, do_err_est);
[tgrid_SO_no_cavity, Pe_SO_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, SO_signal, [0,T_final]);

% COS
[tgrid_COS, Pe_COS, nk_COS, ntot_COS, eps_trunc_COS] = multimode_JC_driven_mode_photon_cap(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, Cos_signal, T_final, do_err_est);
[tgrid_COS_no_cavity, Pe_COS_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, Cos_signal, [0,T_final]);

%% Numerical Errors
fprintf('\n===== Numerical Errors =====\n');
fprintf('SO error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_SO);
fprintf('COS error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_COS);
fprintf('============================\n');

%% PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
plot(tgrid_SO, Pe_SO, 'LineWidth', 5, 'Color', 'red', 'DisplayName', '\textbf{SO}');
plot(tgrid_SO_no_cavity, Pe_SO_no_cavity, '--', 'LineWidth', 3, 'Color', 'k', 'DisplayName', '\textbf{SO no fluc}');
plot(tgrid_COS, Pe_COS, 'LineWidth', 5, 'Color', 'blue', 'DisplayName', '\boldmath$\mathrm{0.9\omega_0}$');
plot(tgrid_COS_no_cavity, Pe_COS_no_cavity, '-.', 'LineWidth', 3, 'Color', 'k', 'DisplayName', '{\boldmath $0.9\omega_0$} \textbf{no fluc}');
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Excitation \ probability}$');
legend('show');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 2: SO Photon Dynamics ---
figure;
hold on;
plot(tgrid_SO, nk_SO, 'LineWidth', 1.5); 
plot(tgrid_SO, ntot_SO, 'k', 'LineWidth', 4, 'DisplayName', '\textbf{Total}'); 
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Photon \ Number}$');
title('\boldmath$\mathrm{SO \ Signal: \ Mode \ Photons \ \& \ Total}$');

% Dynamically generate legend for SO modes in Bold LaTeX
numModesSO = size(nk_SO, 2);
legSO = arrayfun(@(x) sprintf('\\textbf{Mode %d}', x), 1:numModesSO, 'UniformOutput', false);
legSO{end+1} = '\textbf{Total}';
legend(legSO);

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 3: COS Photon Dynamics ---
figure;
hold on;
plot(tgrid_COS, nk_COS, 'LineWidth', 1.5); 
plot(tgrid_COS, ntot_COS, 'k', 'LineWidth', 4, 'DisplayName', '\textbf{Total}'); 
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Photon \ Number}$');
title('\boldmath$\mathrm{COS \ Signal: \ Mode \ Photons \ \& \ Total}$');

% Dynamically generate legend for COS modes in Bold LaTeX
numModesCOS = size(nk_COS, 2);
legCOS = arrayfun(@(x) sprintf('\\textbf{Mode %d}', x), 1:numModesCOS, 'UniformOutput', false);
legCOS{end+1} = '\textbf{Total}';
legend(legCOS);

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);

plot_signals(SO_signal, Cos_signal, angular_freqs_SO,angular_freqs_COS)