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

%% 1. Define Signals
signal_scaling = 7;

[SO_signal, angular_freqs_SO, amps_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, false);
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, false);

% Random phase signals (keeping Baranov amplitudes but with random phase, real-valued)
[Rand_signal_1, ~] = generate_rand_phase(angular_freqs_SO, amps_SO, true, false, 1);
[Rand_signal_2, ~] = generate_rand_phase(angular_freqs_SO, amps_SO, true, false, 2);

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {});

% --- Signal 1: SO ---
sig_configs(1).name = '\textbf{SO}';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;

% --- Signal 2: COS ---
sig_configs(2).name = '\boldmath$\mathbf{0.9\omega_0}$';
sig_configs(2).data = Cos_signal;
sig_configs(2).color = 'b';
sig_configs(2).freqs = angular_freqs_COS;

% --- Signal 3: Rand Phase 1 ---
sig_configs(3).name = '\textbf{Rand Phase 1}';
sig_configs(3).data = Rand_signal_1;
sig_configs(3).color = 'm';
sig_configs(3).freqs = angular_freqs_SO;

% --- Signal 4: Rand Phase 2 ---
sig_configs(4).name = '\textbf{Rand Phase 2}';
sig_configs(4).data = Rand_signal_2;
sig_configs(4).color = [0, 0.5, 0];
sig_configs(4).freqs = angular_freqs_SO;

% Normalize all signals
sig_configs = normalize_sig_configs(sig_configs, 'energy');

% Apply defaults (auto-colors)
sig_configs = prepare_signal_config(sig_configs);


%% 2. Dynamics computation
dt = 0.01;
t_grid = (0:dt:T_final)';

for i = 1:length(sig_configs)
    [tgrid_out, rho_out] = optical_bloch(t_grid, rho_init, nu0, params, sig_configs(i).data);
    sig_configs(i).tgrid = tgrid_out;
    sig_configs(i).Pe = (1 + rho_out(:, 3)) / 2;
end

%% 3. PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
for i = 1:length(sig_configs)
    plot(sig_configs(i).tgrid, sig_configs(i).Pe, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end
grid on;

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Excitation \ probability}$');
legend('show');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);


% --- Figure 2: Full Simulation Instantaneous Frequency & Excitation Analysis ---
% Time grid expanded to cover the full duration: [0, T_final]
t_inst = linspace(0, T_final, 30000)'; 
dt_inst = t_inst(2) - t_inst(1);

for i = 1:length(sig_configs)
    y_sig = arrayfun(sig_configs(i).data, t_inst);
    inst_freq = compute_instantaneous_frequency(y_sig, dt_inst);
    Pe_inst = interp1(sig_configs(i).tgrid, sig_configs(i).Pe, t_inst);
    
    figure('Color', 'w', 'Name', sprintf('%s: Full Analysis', sig_configs(i).name));
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    
    ax(1) = nexttile; 
    plot(t_inst, Pe_inst, 'Color', sig_configs(i).color);
    ylabel('P_e'); title(sprintf('%s: Excitation & Instantaneous Frequency (Full Simulation)', sig_configs(i).name));
    grid on; PlotUtils.styleAxes(gca);
    
    ax(2) = nexttile; 
    hold on;
    plot(t_inst, real(y_sig), 'LineWidth', 1.5, 'Color', sig_configs(i).color, 'DisplayName', 'wave');
    plot(t_inst, inst_freq, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
    xlabel('time'); grid on; legend('Location', 'northeast');
    PlotUtils.styleAxes(gca);
    
    linkaxes(ax, 'x'); 
    xlim([0, T_final]);
end

% Call general plotting utility for the signals
plot_signals(sig_configs);