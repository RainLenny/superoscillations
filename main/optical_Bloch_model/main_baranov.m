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
rho_init = [0; 0; -1]; % Initial conditions of the TLS

params.T1         = 500;
params.T2         = 300;
params.rho30      = -1.0; % Value to which the system relaxes to
params.OmegaTilde = 0.01; % Field-TLS coupling 

T_final = 600;

%% SIGNALS:
signal_scaling = 7;

[SO_signal, angular_freqs_SO, amps_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, false);

[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, false);

% Random phase signals (keeping Baranov amplitudes but with random phase, real-valued)
[Rand_signal_1, ~] = generate_rand_phase(angular_freqs_SO, amps_SO, true, false, 1);
[Rand_signal_2, ~] = generate_rand_phase(angular_freqs_SO, amps_SO, true, false, 2);

% Normalize all signals
[SO_signal, Cos_signal, Rand_signal_1, Rand_signal_2] = normalize_signals({SO_signal, Cos_signal, Rand_signal_1, Rand_signal_2}, 'energy');

%% Dynamics computation
dt = 0.01;
t_grid = (0:dt:T_final)';

% SO
[tgrid_SO, rho_SO] = optical_bloch(t_grid, rho_init, nu0, params, SO_signal);
Pe_SO = (1 + rho_SO(:, 3)) / 2;

% COS
[tgrid_COS, rho_COS] = optical_bloch(t_grid, rho_init, nu0, params, Cos_signal);
Pe_COS = (1 + rho_COS(:, 3)) / 2;

% RAND 1
[tgrid_RAND1, rho_RAND1] = optical_bloch(t_grid, rho_init, nu0, params, Rand_signal_1);
Pe_RAND1 = (1 + rho_RAND1(:, 3)) / 2;

% RAND 2
[tgrid_RAND2, rho_RAND2] = optical_bloch(t_grid, rho_init, nu0, params, Rand_signal_2);
Pe_RAND2 = (1 + rho_RAND2(:, 3)) / 2;


%% PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
plot(tgrid_SO, Pe_SO, 'LineWidth', 5, 'Color', 'red', 'DisplayName', '\textbf{SO}');
plot(tgrid_COS, Pe_COS, 'LineWidth', 5, 'Color', 'blue', 'DisplayName', '\boldmath$\mathrm{0.9\omega_0}$');
plot(tgrid_RAND1, Pe_RAND1, 'LineWidth', 5, 'Color', 'm', 'DisplayName', '\textbf{Rand Phase 1}');
plot(tgrid_RAND2, Pe_RAND2, 'LineWidth', 5, 'Color', [0, 0.5, 0], 'DisplayName', '\textbf{Rand Phase 2}');
grid on;

xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Excitation \ probability}$');
legend('show');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 2: Full Simulation Instantaneous Frequency & Excitation Analysis ---
% Time grid expanded to cover the full duration: [0, T_final]
t_inst = linspace(0, T_final, 30000)'; 
dt_inst = t_inst(2) - t_inst(1);

y_so    = arrayfun(SO_signal, t_inst);
y_cos   = arrayfun(Cos_signal, t_inst);
y_rand1 = arrayfun(Rand_signal_1, t_inst);
y_rand2 = arrayfun(Rand_signal_2, t_inst);

% Compute instantaneous frequency for each signal using the refactored function
inst_freq_so    = compute_instantaneous_frequency(y_so, dt_inst);
inst_freq_cos   = compute_instantaneous_frequency(y_cos, dt_inst);
inst_freq_rand1 = compute_instantaneous_frequency(y_rand1, dt_inst);
inst_freq_rand2 = compute_instantaneous_frequency(y_rand2, dt_inst);

% Interpolate Excitation Probabilities to match the expanded t_inst grid precisely
Pe_SO_inst    = interp1(tgrid_SO, Pe_SO, t_inst);
Pe_COS_inst   = interp1(tgrid_COS, Pe_COS, t_inst);
Pe_RAND1_inst = interp1(tgrid_RAND1, Pe_RAND1, t_inst);
Pe_RAND2_inst = interp1(tgrid_RAND2, Pe_RAND2, t_inst);


% --- Figure 2.1: SO ---
figure('Color', 'w', 'Name', 'SO: Full Analysis');
tiledlayout(2, 1, 'TileSpacing', 'compact');

ax_so(1) = nexttile; 
plot(t_inst, Pe_SO_inst, 'LineWidth', 2.5, 'Color', 'r');
ylabel('P_e'); title('SO: Excitation & Instantaneous Frequency (Full Simulation)');
grid on; PlotUtils.styleAxes(gca);

ax_so(2) = nexttile; 
hold on;
plot(t_inst, real(y_so), 'LineWidth', 1.5, 'Color', 'r', 'DisplayName', 'wave');
plot(t_inst, inst_freq_so, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
xlabel('time'); grid on; legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

linkaxes(ax_so, 'x'); 
xlim([0, T_final]);


% --- Figure 2.2: COS ---
figure('Color', 'w', 'Name', 'COS: Full Analysis');
tiledlayout(2, 1, 'TileSpacing', 'compact');

ax_cos(1) = nexttile; 
plot(t_inst, Pe_COS_inst, 'LineWidth', 2.5, 'Color', 'b');
ylabel('P_e'); title('COS (0.9): Excitation & Instantaneous Frequency (Full Simulation)');
grid on; PlotUtils.styleAxes(gca);

ax_cos(2) = nexttile; 
hold on;
plot(t_inst, real(y_cos), 'LineWidth', 1.5, 'Color', 'b', 'DisplayName', 'wave');
plot(t_inst, inst_freq_cos, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
xlabel('time'); grid on; legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

linkaxes(ax_cos, 'x'); 
xlim([0, T_final]);


% --- Figure 2.3: RAND PHASE 1 ---
figure('Color', 'w', 'Name', 'Rand Phase 1: Full Analysis');
tiledlayout(2, 1, 'TileSpacing', 'compact');

ax_rand1(1) = nexttile; 
plot(t_inst, Pe_RAND1_inst, 'LineWidth', 2.5, 'Color', 'm');
ylabel('P_e'); title('Rand Phase 1: Excitation & Instantaneous Frequency (Full Simulation)');
grid on; PlotUtils.styleAxes(gca);

ax_rand1(2) = nexttile; 
hold on;
plot(t_inst, real(y_rand1), 'LineWidth', 1.5, 'Color', 'm', 'DisplayName', 'wave');
plot(t_inst, inst_freq_rand1, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
xlabel('time'); grid on; legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

linkaxes(ax_rand1, 'x'); 
xlim([0, T_final]);


% --- Figure 2.4: RAND PHASE 2 ---
figure('Color', 'w', 'Name', 'Rand Phase 2: Full Analysis');
tiledlayout(2, 1, 'TileSpacing', 'compact');

ax_rand2(1) = nexttile; 
plot(t_inst, Pe_RAND2_inst, 'LineWidth', 2.5, 'Color', [0, 0.5, 0]);
ylabel('P_e'); title('Rand Phase 2: Excitation & Instantaneous Frequency (Full Simulation)');
grid on; PlotUtils.styleAxes(gca);

ax_rand2(2) = nexttile; 
hold on;
plot(t_inst, real(y_rand2), 'LineWidth', 1.5, 'Color', [0, 0.5, 0], 'DisplayName', 'wave');
plot(t_inst, inst_freq_rand2, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
xlabel('time'); grid on; legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

linkaxes(ax_rand2, 'x'); 
xlim([0, T_final]);


% Call general plotting utility for the four signals
plot_signals(SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS, Rand_signal_1, Rand_signal_2);