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
[Rand_signal, ~] = generate_rand_phase(angular_freqs_SO, amps_baranov, true, false, 5);

% Normalize all signals symbolically by the peak of the first signal
[SO_signal, Cos_signal, Flat_signal, Rand_signal] = normalize_signals({SO_signal, Cos_signal, Flat_signal, Rand_signal}, 'peak');

%% Dynamics computation
dt = 0.01;
t_grid = (0:dt:T_final)';

% SO
[tgrid_SO, rho_SO] = optical_bloch(t_grid, rho_init, nu0, params, SO_signal);
Pe_SO = (1 + rho_SO(:, 3)) / 2;

% COS
[tgrid_COS, rho_COS] = optical_bloch(t_grid, rho_init, nu0, params, Cos_signal);
Pe_COS = (1 + rho_COS(:, 3)) / 2;

% FLAT
[tgrid_FLAT, rho_FLAT] = optical_bloch(t_grid, rho_init, nu0, params, Flat_signal);
Pe_FLAT = (1 + rho_FLAT(:, 3)) / 2;

% RAND
[tgrid_RAND, rho_RAND] = optical_bloch(t_grid, rho_init, nu0, params, Rand_signal);
Pe_RAND = (1 + rho_RAND(:, 3)) / 2;


%% PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
plot(tgrid_SO, Pe_SO, 'LineWidth', 5, 'Color', 'red', 'DisplayName', '\textbf{SO}');
plot(tgrid_COS, Pe_COS, 'LineWidth', 5, 'Color', 'blue', 'DisplayName', '\boldmath$\mathrm{0.9\omega_0}$');
plot(tgrid_FLAT, Pe_FLAT, 'LineWidth', 5, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Flat}');
plot(tgrid_RAND, Pe_RAND, 'LineWidth', 5, 'Color', 'm', 'DisplayName', '\textbf{Rand Phase}');
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Excitation \ probability}$');
legend('show');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 2: Instantaneous Frequency Analysis ---
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

% --- Figure 1: SO ---
figure('Color', 'w', 'Name', 'SO: Instantaneous Frequency Analysis');
hold on;
plot(t_inst, real(y_so), 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'wave');
plot(t_inst, inst_freq_so, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('SO: Wave and Instantaneous Frequency');
xlabel('time');
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% --- Figure 2: COS ---
figure('Color', 'w', 'Name', 'COS: Instantaneous Frequency Analysis');
hold on;
plot(t_inst, real(y_cos), 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'wave');
plot(t_inst, inst_freq_cos, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('COS (0.9): Wave and Instantaneous Frequency');
xlabel('time');
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% --- Figure 3: FLAT ---
figure('Color', 'w', 'Name', 'Flat Spectrum: Instantaneous Frequency Analysis');
hold on;
plot(t_inst, real(y_flat), 'LineWidth', 2, 'Color', [0, 0.5, 0], 'DisplayName', 'wave');
plot(t_inst, inst_freq_flat, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('Flat Spectrum: Wave and Instantaneous Frequency');
xlabel('time');
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% --- Figure 4: RAND PHASE ---
figure('Color', 'w', 'Name', 'Rand Phase: Instantaneous Frequency Analysis');
hold on;
plot(t_inst, real(y_rand), 'LineWidth', 2, 'Color', 'm', 'DisplayName', 'wave');
plot(t_inst, inst_freq_rand, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('Rand Phase: Wave and Instantaneous Frequency');
xlabel('time');
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% Call general plotting utility for the four signals
plot_signals(SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS, Flat_signal, Rand_signal);