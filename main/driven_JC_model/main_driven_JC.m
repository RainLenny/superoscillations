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
J_drive = 3;
J_fluc  = 0.006;

T_final = 500;
nmax    = 5;

% Control flags
do_err_est = 0;

%% 1. Define Signals
signal_scaling = 1.3;

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {});

% --- Signal 1: SO ---
[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, true);
sig_configs(1).name = '\textbf{SO }';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;

% --- Signal 2: COS ---
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, true);
sig_configs(2).name = '\boldmath$\mathbf{0.9\omega_0}$';
sig_configs(2).data = Cos_signal;
sig_configs(2).color = 'b';
sig_configs(2).freqs = angular_freqs_COS;

% % --- Signal 3: FLAT ---
% [Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, true);
% sig_configs(3).name = '\textbf{Flat}';
% sig_configs(3).data = Flat_signal;
% sig_configs(3).color = [0, 0.5, 0];
% sig_configs(3).freqs = angular_freqs_SO;
% 
% --- Signal 4: RAND ---
% data_baranov = load('SO_Baranov.mat', 'amps_SO');
% amps_baranov = data_baranov.amps_SO * signal_scaling;
% [Rand_signal, ~] = generate_rand_phase(angular_freqs_SO, amps_baranov, true, true, 4);
% sig_configs(3).name = '\textbf{Rand Phase}';
% sig_configs(3).data = Rand_signal;
% sig_configs(3).color = 'g';
% sig_configs(3).freqs = angular_freqs_SO;

% Normalize all signals symbolically by the peak of the first signal
sigs = {sig_configs.data};
[norm_sigs{1:length(sigs)}] = normalize_signals(sigs, 'energy');
for i = 1:length(sig_configs)
    sig_configs(i).data = norm_sigs{i};
end

% Apply defaults (auto-colors)
sig_configs = prepare_signal_config(sig_configs);

%% 2. Dynamics computation
fprintf('\n===== Numerical Errors =====\n');
for i = 1:length(sig_configs)
    [tgrid, Pe, nk, ntot, eps_trunc] = multimode_JC_driven_mode_photon_cap(sig_configs(i).freqs, nu0, nmax, J_fluc, J_drive, sig_configs(i).data, T_final, do_err_est);
    [tgrid_nc, Pe_nc] = JC_drive_only(nu0, J_fluc*J_drive, sig_configs(i).data, [0, T_final]);
    
    sig_configs(i).tgrid = tgrid;
    sig_configs(i).Pe = Pe;
    sig_configs(i).nk = nk;
    sig_configs(i).ntot = ntot;
    sig_configs(i).eps_trunc = eps_trunc;
    
    sig_configs(i).tgrid_nc = tgrid_nc;
    sig_configs(i).Pe_nc = Pe_nc;
    
    fprintf('%s error (Hilbert truncation nmax+1): %.8e\n', strrep(sig_configs(i).name, '\', ''), eps_trunc);
end
fprintf('============================\n');

%% 3. PLOTS
PlotUtils.setupDefaults();

% --- Figure 1: Original Excitation Probability ---
figure;
hold on;
for i = 1:length(sig_configs)
    plot(sig_configs(i).tgrid, sig_configs(i).Pe, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
    plot(sig_configs(i).tgrid_nc, sig_configs(i).Pe_nc, '--', 'Color', 'k', 'DisplayName', [sig_configs(i).name, ' \textbf{no fluc}']);
end
grid on;


xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Excitation \ probability}$');
legend('show');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);

% --- Figure 1b: Excitation Probability (no fluc lines removed) ---
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

% --- Figure 2 to N+1: Photon Dynamics for each signal ---
for i = 1:length(sig_configs)
    figure;
    hold on;
    plot(sig_configs(i).tgrid, sig_configs(i).nk, 'LineWidth', 1.5); 
    plot(sig_configs(i).tgrid, sig_configs(i).ntot, 'k', 'DisplayName', '\textbf{Total}'); 
    grid on;

    xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
    ylabel('\boldmath$\mathbf{Photon \ Number}$');
    title(sprintf('%s \\textbf{ Signal: Mode Photons \\& Total}', sig_configs(i).name));

    % Dynamically generate legend for modes in Bold LaTeX
    numModes = size(sig_configs(i).nk, 2);
    leg = arrayfun(@(x) sprintf('\\textbf{Mode %d}', x), 1:numModes, 'UniformOutput', false);
    leg{end+1} = '\textbf{Total}';
    legend(leg);

    % Apply Bold LaTeX Axis Ticks
    PlotUtils.styleAxes(gca);
end


% --- Figure N+2: Full Simulation Instantaneous Frequency & Excitation Analysis ---
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