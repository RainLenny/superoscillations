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

do_err_est = 0;

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

%% PLOT

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
plot(tgrid_SO, Pe_SO, 'LineWidth', 5, 'Color', 'red');
plot(tgrid_SO_no_cavity, Pe_SO_no_cavity, '--', 'LineWidth', 3, 'Color', 'k');
plot(tgrid_COS, Pe_COS, 'LineWidth', 5, 'Color', 'blue');
plot(tgrid_COS_no_cavity, Pe_COS_no_cavity, '-.', 'LineWidth', 3, 'Color', 'k');
grid on;
xlabel('Time [arb]');
ylabel('Excitation probability');
ax1 = gca; % Get current axes
ax1.FontWeight = 'bold';
ax1.FontSize   = 12;
legend({'SO','SO no fluc','COS','COS no fluc'}, 'Location','best');

% --- Figure 2: SO Photon Dynamics ---
figure;
hold on;
plot(tgrid_SO, nk_SO, 'LineWidth', 1.5); % Plots each mode
plot(tgrid_SO, ntot_SO, 'k', 'LineWidth', 3); % Overlay Total Photons in black
grid on;
xlabel('Time [arb]');
ylabel('Photon Number');
title('SO Signal: Mode Photons & Total');
ax2 = gca; % Get current axes
ax2.FontWeight = 'bold';
ax2.FontSize   = 12;

% Dynamically generate legend for SO modes + Total
numModesSO = size(nk_SO, 2); % Assumes modes are in columns
legSO = arrayfun(@(x) sprintf('Mode %d', x), 1:numModesSO, 'UniformOutput', false);
legSO{end+1} = 'Total';
legend(legSO, 'Location', 'best');


% --- Figure 3: COS Photon Dynamics ---
figure;
hold on;
plot(tgrid_COS, nk_COS, 'LineWidth', 1.5); % Plots each mode
plot(tgrid_COS, ntot_COS, 'k', 'LineWidth', 3); % Overlay Total Photons in black
grid on;
xlabel('Time [arb]');
ylabel('Photon Number');
title('COS Signal: Mode Photons & Total');
ax3 = gca; % Get current axes
ax3.FontWeight = 'bold';
ax3.FontSize   = 12;

% Dynamically generate legend for COS modes + Total
numModesCOS = size(nk_COS, 2); % Assumes modes are in columns
legCOS = arrayfun(@(x) sprintf('Mode %d', x), 1:numModesCOS, 'UniformOutput', false);
legCOS{end+1} = 'Total';
legend(legCOS, 'Location', 'best');

% plot_signals(freq_scaling,amp_scaling)