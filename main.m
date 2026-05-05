%% CONSTANTS
clear; clc
freq_scaling = 1;
amp_scaling = 1;
[SO_signal, Cos_signal,angular_freqs_SO, angular_freqs_COS]  = generate_signals_Baranov_2014(freq_scaling,amp_scaling);
nu0 = 1;
J_drive = 6;
J_fluc = 0.01;

T_final = 500;
nmax = 3;

do_err_est = 1;

REMOVE MINUS SIGN FROM EFFECTIVE DRIVE - MATCH TO PAPER 

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

%% PLOT STYLED WITH BOLD LATEX

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
plot(tgrid_SO, Pe_SO, 'LineWidth', 5, 'Color', 'red', 'DisplayName', '\textbf{SO}');
plot(tgrid_SO_no_cavity, Pe_SO_no_cavity, '--', 'LineWidth', 3, 'Color', 'k', 'DisplayName', '\textbf{SO no fluc}');
plot(tgrid_COS, Pe_COS, 'LineWidth', 5, 'Color', 'blue', 'DisplayName', '\boldmath$\mathrm{0.9\nu_0}$');
plot(tgrid_COS_no_cavity, Pe_COS_no_cavity, '-.', 'LineWidth', 3, 'Color', 'k', 'DisplayName', '{\boldmath $0.9\nu_0$} \textbf{no fluc}');
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\nu_0]}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('\boldmath$\mathrm{Excitation \ probability}$', 'FontSize', 14, 'Interpreter', 'latex');
legend('show', 'FontSize', 12, 'Location', 'best', 'Interpreter', 'latex', 'FontWeight', 'bold');

% Apply Bold LaTeX Axis Ticks
ax1 = gca;
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 13;
ax1.FontWeight = 'bold';
ax1.XTickLabel = arrayfun(@(x) sprintf('$\\mathbf{%g}$', x), ax1.XTick, 'UniformOutput', false);
ax1.YTickLabel = arrayfun(@(y) sprintf('$\\mathbf{%g}$', y), ax1.YTick, 'UniformOutput', false);


% --- Figure 2: SO Photon Dynamics ---
figure;
hold on;
plot(tgrid_SO, nk_SO, 'LineWidth', 1.5); 
plot(tgrid_SO, ntot_SO, 'k', 'LineWidth', 4, 'DisplayName', '\textbf{Total}'); 
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\nu_0]}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('\boldmath$\mathrm{Photon \ Number}$', 'FontSize', 14, 'Interpreter', 'latex');
title('\boldmath$\mathrm{SO \ Signal: \ Mode \ Photons \ \& \ Total}$', 'Interpreter', 'latex', 'FontSize', 14);

% Dynamically generate legend for SO modes in Bold LaTeX
numModesSO = size(nk_SO, 2);
legSO = arrayfun(@(x) sprintf('\\textbf{Mode %d}', x), 1:numModesSO, 'UniformOutput', false);
legSO{end+1} = '\textbf{Total}';
legend(legSO, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

% Apply Bold LaTeX Axis Ticks
ax2 = gca;
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = 13;
ax2.FontWeight = 'bold';
ax2.XTickLabel = arrayfun(@(x) sprintf('$\\mathbf{%g}$', x), ax2.XTick, 'UniformOutput', false);
ax2.YTickLabel = arrayfun(@(y) sprintf('$\\mathbf{%g}$', y), ax2.YTick, 'UniformOutput', false);


% --- Figure 3: COS Photon Dynamics ---
figure;
hold on;
plot(tgrid_COS, nk_COS, 'LineWidth', 1.5); 
plot(tgrid_COS, ntot_COS, 'k', 'LineWidth', 4, 'DisplayName', '\textbf{Total}'); 
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\nu_0]}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('\boldmath$\mathrm{Photon \ Number}$', 'FontSize', 14, 'Interpreter', 'latex');
title('\boldmath$\mathrm{COS \ Signal: \ Mode \ Photons \ \& \ Total}$', 'Interpreter', 'latex', 'FontSize', 14);

% Dynamically generate legend for COS modes in Bold LaTeX
numModesCOS = size(nk_COS, 2);
legCOS = arrayfun(@(x) sprintf('\\textbf{Mode %d}', x), 1:numModesCOS, 'UniformOutput', false);
legCOS{end+1} = '\textbf{Total}';
legend(legCOS, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

% Apply Bold LaTeX Axis Ticks
ax3 = gca;
ax3.TickLabelInterpreter = 'latex';
ax3.FontSize = 13;
ax3.FontWeight = 'bold';
ax3.XTickLabel = arrayfun(@(x) sprintf('$\\mathbf{%g}$', x), ax3.XTick, 'UniformOutput', false);
ax3.YTickLabel = arrayfun(@(y) sprintf('$\\mathbf{%g}$', y), ax3.YTick, 'UniformOutput', false);

% plot_signals(freq_scaling,amp_scaling)