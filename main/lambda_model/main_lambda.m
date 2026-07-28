%% PATH SETUP & CONSTANTS
clear; clc;

% NEED TO CALIBRATE THE LAMBDA SYSTEM FOR THE FREQUENCIES OF THE SIGNAL!!!!!!!!!
%% CONSTANTS
T_final = 500;
t_0 = 250;
T = 100;

%% 1. Define Signals
signal_scaling = 1;
sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {});

% --- Signal 1: SO ---
[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Derek', 1, signal_scaling, true, false);
sig_configs(1).name = '\textbf{SO}';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;

% --- Signal 2: Flat ---
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, false);
sig_configs(2).name = '\textbf{Flat}';
sig_configs(2).data = Flat_signal;
sig_configs(2).color = [0, 0.5, 0];
sig_configs(2).freqs = angular_freqs_SO;

% Normalize signals together
sigs = {sig_configs.data};
[norm_sigs{1:length(sigs)}] = normalize_signals(sigs, 'peak');
for i = 1:length(sig_configs)
    sig_configs(i).data = norm_sigs{i};
end

% Apply defaults (auto-colors)
sig_configs = prepare_signal_config(sig_configs);

% Define Probe fields Omega_p
Omega_p_max = 5; % Peak amplitude
for i = 1:length(sig_configs)
    sig_configs(i).Omega_p = @(t) Omega_p_max * sig_configs(i).data(t);
end

% Define Stokes field Omega_s (Adiabatic Anchor)
Omega_s0 = 10; % Strong Stokes field
Omega_s_func = @(t) Omega_s0 * exp(-(t - t_0).^2 / (2.5*T)^2);

%% 2. DYNAMICS
delta = 0; % Two-photon resonance
Delta = 0; % Single-photon detuning
c_init = [1; 0; 0]; % Start in |1>
t_span = [0, T_final];

for i = 1:length(sig_configs)
    [tgrid, rho] = lambda_system_dynamics(t_span, c_init, sig_configs(i).Omega_p, Omega_s_func, delta, Delta);
    sig_configs(i).tgrid = tgrid;
    sig_configs(i).rho = rho;
end

%% 3. PLOTS
PlotUtils.setupDefaults();

% =========================================================================
% DYNAMICS PLOTS FOR EACH SIGNAL
% =========================================================================

for i = 1:length(sig_configs)
    % --- Driving Fields and Populations ---
    figure('Color', 'w', 'Name', sprintf('Lambda System Dynamics (%s Signal)', sig_configs(i).name));

    subplot(2,1,1);
    hold on;
    plot(sig_configs(i).tgrid, abs(sig_configs(i).Omega_p(sig_configs(i).tgrid)), 'Color', sig_configs(i).color, 'LineWidth', 2, 'DisplayName', sprintf('|\\Omega_p(t)| (Probe - %s)', sig_configs(i).name));
    plot(sig_configs(i).tgrid, abs(Omega_s_func(sig_configs(i).tgrid)), 'b', 'LineWidth', 2, 'DisplayName', '|\Omega_s(t)| (Stokes - Anchor)');
    grid on;
    xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
    ylabel('\boldmath$\mathrm{Rabi \ Frequencies}$');
    title(sprintf('\\boldmath$\\mathrm{Driving \\ Fields \\ (%s \\ Signal)}$', sig_configs(i).name));
    legend('show');
    PlotUtils.styleAxes(gca);

    subplot(2,1,2);
    hold on;
    plot(sig_configs(i).tgrid, sig_configs(i).rho(:,1), 'LineWidth', 2, 'DisplayName', '\boldmath$\rho_{11} \ \mathrm{(Ground \ 1)}$');
    plot(sig_configs(i).tgrid, sig_configs(i).rho(:,2), 'LineWidth', 2, 'DisplayName', '\boldmath$\rho_{22} \ \mathrm{(Ground \ 2)}$');
    plot(sig_configs(i).tgrid, sig_configs(i).rho(:,3), 'k', 'DisplayName', '\boldmath$\rho_{33} \ \mathrm{(Excited)}$');
    grid on;
    xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
    ylabel('\boldmath$\mathrm{Populations}$');
    title(sprintf('\\boldmath$\\mathrm{State \\ Populations \\ (%s \\ Signal)}$', sig_configs(i).name));
    legend('show');
    PlotUtils.styleAxes(gca);

    % --- Zoom on Excited State Population ---
    figure('Color', 'w', 'Name', sprintf('Diabatic Transition Zoom (%s Signal)', sig_configs(i).name));
    hold on;
    plot(sig_configs(i).tgrid, sig_configs(i).rho(:,3), 'k', 'DisplayName', '\boldmath$\rho_{33} \ \mathrm{(Excited)}$');
    grid on;
    xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
    ylabel('\boldmath$\mathrm{Population \ \rho_{33}}$');
    title(sprintf('\\boldmath$\\mathrm{Diabatic \\ Transition \\ to \\ Excited \\ State \\ (%s \\ Signal)}$', sig_configs(i).name));
    xlim([245, 255]); % Zoom in precisely on the region
    legend('show');
    PlotUtils.styleAxes(gca);
end

% =========================================================================
% SIGNAL ANALYSIS: TIME DOMAIN, FFT, AND INSTANTANEOUS FREQUENCY
% =========================================================================

% We can just call plot_signals which does Time Domain and FFT!
plot_signals(sig_configs);

% --- Instantaneous Frequency Analysis ---
t_inst = linspace(200, 300, 5000)';
dt_inst = t_inst(2) - t_inst(1);

figure('Color', 'w', 'Name', 'Instantaneous Frequency Analysis');

x_limits = [245, 255];
y_limits = [-2, 4];
num_sigs = length(sig_configs);

for i = 1:num_sigs
    y_inst = arrayfun(sig_configs(i).data, t_inst);
    inst_freq = compute_instantaneous_frequency(y_inst, dt_inst);

    subplot(num_sigs, 1, i);
    hold on;
    plot(t_inst, real(y_inst), 'LineWidth', 2, 'Color', sig_configs(i).color, 'DisplayName', 'Wave');
    plot(t_inst, inst_freq, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
    title(sprintf('%s: Wave and Instantaneous Frequency', sig_configs(i).name));
    xlabel('time');
    xlim(x_limits);
    ylim(y_limits);
    grid on;
    legend('Location', 'northeast');
    PlotUtils.styleAxes(gca);
end
