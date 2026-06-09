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

[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, true);

[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9,signal_scaling,true,true);

% Flat spectrum signal with same frequencies
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, true);

% Random phase signal (like the one in main_filters_noise.m, keeping Baranov amplitudes)
N_SO = length(angular_freqs_SO);
data_baranov = load('SO_Baranov.mat', 'amps_SO');
amps_baranov = data_baranov.amps_SO * signal_scaling;
rng(1);
tau_rand = rand(1, N_SO);
amps_rand = abs(amps_baranov) .* exp(-1i .* tau_rand);
[Rand_signal, ~] = generate_signal_base(angular_freqs_SO, amps_rand, true, true);

% Normalize all signals symbolically by the peak of the first signal
[SO_signal, Cos_signal, Flat_signal, Rand_signal] = normalize_signals({SO_signal, Cos_signal, Flat_signal, Rand_signal}, 'peak');

%% Dynamics computation
% SO
[tgrid_SO, Pe_SO, nk_SO, ntot_SO, eps_trunc_SO] = multimode_JC_driven_mode_photon_cap(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, SO_signal, T_final, do_err_est);
[tgrid_SO_no_cavity, Pe_SO_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, SO_signal, [0,T_final]);

% COS
[tgrid_COS, Pe_COS, nk_COS, ntot_COS, eps_trunc_COS] = multimode_JC_driven_mode_photon_cap(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, Cos_signal, T_final, do_err_est);
[tgrid_COS_no_cavity, Pe_COS_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, Cos_signal, [0,T_final]);

% FLAT
[tgrid_FLAT, Pe_FLAT, nk_FLAT, ntot_FLAT, eps_trunc_FLAT] = multimode_JC_driven_mode_photon_cap(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, Flat_signal, T_final, do_err_est);
[tgrid_FLAT_no_cavity, Pe_FLAT_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, Flat_signal, [0,T_final]);

% RAND
[tgrid_RAND, Pe_RAND, nk_RAND, ntot_RAND, eps_trunc_RAND] = multimode_JC_driven_mode_photon_cap(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, Rand_signal, T_final, do_err_est);
[tgrid_RAND_no_cavity, Pe_RAND_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, Rand_signal, [0,T_final]);

%% Numerical Errors
fprintf('\n===== Numerical Errors =====\n');
fprintf('SO error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_SO);
fprintf('COS error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_COS);
fprintf('FLAT error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_FLAT);
fprintf('RAND error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_RAND);
fprintf('============================\n');

%% PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
plot(tgrid_SO, Pe_SO, 'LineWidth', 5, 'Color', 'red', 'DisplayName', '\textbf{SO}');
plot(tgrid_SO_no_cavity, Pe_SO_no_cavity, '--', 'LineWidth', 3, 'Color', 'red', 'DisplayName', '\textbf{SO no fluc}');
plot(tgrid_COS, Pe_COS, 'LineWidth', 5, 'Color', 'blue', 'DisplayName', '\boldmath$\mathrm{0.9\omega_0}$');
plot(tgrid_COS_no_cavity, Pe_COS_no_cavity, '-.', 'LineWidth', 3, 'Color', 'blue', 'DisplayName', '{\boldmath $0.9\omega_0$} \textbf{no fluc}');
plot(tgrid_FLAT, Pe_FLAT, 'LineWidth', 5, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');
plot(tgrid_FLAT_no_cavity, Pe_FLAT_no_cavity, ':', 'LineWidth', 3, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat no fluc}');
plot(tgrid_RAND, Pe_RAND, 'LineWidth', 5, 'Color', 'm', 'DisplayName', '\textbf{Rand Phase}');
plot(tgrid_RAND_no_cavity, Pe_RAND_no_cavity, '--', 'LineWidth', 3, 'Color', 'm', 'DisplayName', '\textbf{Rand Phase no fluc}');
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


% --- Figure 4: FLAT Photon Dynamics ---
figure;
hold on;
plot(tgrid_FLAT, nk_FLAT, 'LineWidth', 1.5); 
plot(tgrid_FLAT, ntot_FLAT, 'k', 'LineWidth', 4, 'DisplayName', '\textbf{Total}'); 
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Photon \ Number}$');
title('\boldmath$\mathrm{FLAT \ Signal: \ Mode \ Photons \ \& \ Total}$');

% Dynamically generate legend for FLAT modes in Bold LaTeX
numModesFLAT = size(nk_FLAT, 2);
legFLAT = arrayfun(@(x) sprintf('\\textbf{Mode %d}', x), 1:numModesFLAT, 'UniformOutput', false);
legFLAT{end+1} = '\textbf{Total}';
legend(legFLAT);

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 5: RAND Photon Dynamics ---
figure;
hold on;
plot(tgrid_RAND, nk_RAND, 'LineWidth', 1.5); 
plot(tgrid_RAND, ntot_RAND, 'k', 'LineWidth', 4, 'DisplayName', '\textbf{Total}'); 
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Photon \ Number}$');
title('\boldmath$\mathrm{RAND \ Signal: \ Mode \ Photons \ \& \ Total}$');

% Dynamically generate legend for RAND modes in Bold LaTeX
numModesRAND = size(nk_RAND, 2);
legRAND = arrayfun(@(x) sprintf('\\textbf{Mode %d}', x), 1:numModesRAND, 'UniformOutput', false);
legRAND{end+1} = '\textbf{Total}';
legend(legRAND);

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 6: Instantaneous Frequency Analysis ---
% Time grid for instantaneous frequency analysis
t_inst = linspace(200, 300, 5000)';
dt_inst = t_inst(2) - t_inst(1);

y_so   = arrayfun(SO_signal, t_inst);
y_cos  = arrayfun(Cos_signal, t_inst);
y_flat = arrayfun(Flat_signal, t_inst);
y_rand = arrayfun(Rand_signal, t_inst);

% Compute instantaneous frequency for each signal using the refactored function
inst_freq_so = compute_instantaneous_frequency(y_so, dt_inst);
inst_freq_cos = compute_instantaneous_frequency(y_cos, dt_inst);
inst_freq_flat = compute_instantaneous_frequency(y_flat, dt_inst);
inst_freq_rand = compute_instantaneous_frequency(y_rand, dt_inst);

figure('Color', 'w', 'Name', 'Instantaneous Frequency Analysis');

x_limits = [245, 255];
y_limits = [-2, 4];

% Subplot 1: SO
subplot(4,1,1);
hold on;
plot(t_inst, real(y_so), 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'wave');
plot(t_inst, inst_freq_so, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('SO: Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% Subplot 2: COS
subplot(4,1,2);
hold on;
plot(t_inst, real(y_cos), 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'wave');
plot(t_inst, inst_freq_cos, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('COS (0.9): Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% Subplot 3: FLAT
subplot(4,1,3);
hold on;
plot(t_inst, real(y_flat), 'LineWidth', 2, 'Color', [0, 0.5, 0], 'DisplayName', 'wave');
plot(t_inst, inst_freq_flat, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('Flat Spectrum: Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% Subplot 4: RAND PHASE
subplot(4,1,4);
hold on;
plot(t_inst, real(y_rand), 'LineWidth', 2, 'Color', 'm', 'DisplayName', 'wave');
plot(t_inst, inst_freq_rand, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('Rand Phase: Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

plot_signals(SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS, Flat_signal, Rand_signal)