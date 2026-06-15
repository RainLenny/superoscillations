%% PATH SETUP & CONSTANTS
clear; clc;
% Add all project subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

%% CONSTANTS
T_final = 500;
t_0 = 250;
T = 100;

%% SIGNALS
signal_scaling = 1;

% 1. Generate Superoscillatory (SO) signal
[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, true);

% 2. Generate Flat spectrum signal
[Flat_signal, ~] = generate_equal_spread(angular_freqs_SO, 1, signal_scaling, true, true);

% 3. Normalize signals together
[SO_signal, Flat_signal] = normalize_signals({SO_signal, Flat_signal}, 'peak');

% Define Probe fields Omega_p
Omega_p_max = 5; % Peak amplitude
Omega_p_SO = @(t) Omega_p_max * SO_signal(t);
Omega_p_Flat = @(t) Omega_p_max * Flat_signal(t);

% Define Stokes field Omega_s (Adiabatic Anchor)
Omega_s0 = 10; % Strong Stokes field
Omega_s_func = @(t) Omega_s0 * exp(-(t - t_0).^2 / (2.5*T)^2);

%% DYNAMICS
delta = 0; % Two-photon resonance
Delta = 0; % Single-photon detuning
c_init = [1; 0; 0]; % Start in |1>
t_span = [0, T_final];

% Run the Lambda system simulation for SO signal
[tgrid_SO, rho_SO] = lambda_system_dynamics(t_span, c_init, Omega_p_SO, Omega_s_func, delta, Delta);

% Run the Lambda system simulation for Flat signal
[tgrid_Flat, rho_Flat] = lambda_system_dynamics(t_span, c_init, Omega_p_Flat, Omega_s_func, delta, Delta);

%% PLOTS
PlotUtils.setupDefaults();

% =========================================================================
% 1. SO SIGNAL DYNAMICS PLOTS
% =========================================================================

% --- Figure 1: SO Signal Driving Fields and Populations ---
figure('Color', 'w', 'Name', 'Lambda System Dynamics (SO Signal)');

subplot(2,1,1);
hold on;
plot(tgrid_SO, abs(Omega_p_SO(tgrid_SO)), 'r', 'LineWidth', 2, 'DisplayName', '|\Omega_p(t)| (Probe - SO)');
plot(tgrid_SO, abs(Omega_s_func(tgrid_SO)), 'b', 'LineWidth', 2, 'DisplayName', '|\Omega_s(t)| (Stokes - Anchor)');
grid on;
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Rabi \ Frequencies}$');
title('\boldmath$\mathrm{Driving \ Fields \ (SO \ Signal)}$');
legend('show');
PlotUtils.styleAxes(gca);

subplot(2,1,2);
hold on;
plot(tgrid_SO, rho_SO(:,1), 'LineWidth', 2, 'DisplayName', '\boldmath$\rho_{11} \ \mathrm{(Ground \ 1)}$');
plot(tgrid_SO, rho_SO(:,2), 'LineWidth', 2, 'DisplayName', '\boldmath$\rho_{22} \ \mathrm{(Ground \ 2)}$');
plot(tgrid_SO, rho_SO(:,3), 'k', 'LineWidth', 3, 'DisplayName', '\boldmath$\rho_{33} \ \mathrm{(Excited)}$');
grid on;
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Populations}$');
title('\boldmath$\mathrm{State \ Populations \ (SO \ Signal)}$');
legend('show');
PlotUtils.styleAxes(gca);

% --- Figure 2: Zoom on SO Excited State Population ---
figure('Color', 'w', 'Name', 'Diabatic Transition Zoom (SO Signal)');
hold on;
plot(tgrid_SO, rho_SO(:,3), 'k', 'LineWidth', 4, 'DisplayName', '\boldmath$\rho_{33} \ \mathrm{(Excited)}$');
grid on;
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Population \ \rho_{33}}$');
title('\boldmath$\mathrm{Diabatic \ Transition \ to \ Excited \ State \ (SO \ Signal)}$');
xlim([245, 255]); % Zoom in precisely on the superoscillatory region
legend('show');
PlotUtils.styleAxes(gca);

% =========================================================================
% 2. FLAT SIGNAL DYNAMICS PLOTS
% =========================================================================

% --- Figure 3: Flat Signal Driving Fields and Populations ---
figure('Color', 'w', 'Name', 'Lambda System Dynamics (Flat Signal)');

subplot(2,1,1);
hold on;
plot(tgrid_Flat, abs(Omega_p_Flat(tgrid_Flat)), 'Color', [0, 0.5, 0], 'LineWidth', 2, 'DisplayName', '|\Omega_p(t)| (Probe - Flat)');
plot(tgrid_Flat, abs(Omega_s_func(tgrid_Flat)), 'b', 'LineWidth', 2, 'DisplayName', '|\Omega_s(t)| (Stokes - Anchor)');
grid on;
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Rabi \ Frequencies}$');
title('\boldmath$\mathrm{Driving \ Fields \ (Flat \ Signal)}$');
legend('show');
PlotUtils.styleAxes(gca);

subplot(2,1,2);
hold on;
plot(tgrid_Flat, rho_Flat(:,1), 'LineWidth', 2, 'DisplayName', '\boldmath$\rho_{11} \ \mathrm{(Ground \ 1)}$');
plot(tgrid_Flat, rho_Flat(:,2), 'LineWidth', 2, 'DisplayName', '\boldmath$\rho_{22} \ \mathrm{(Ground \ 2)}$');
plot(tgrid_Flat, rho_Flat(:,3), 'k', 'LineWidth', 3, 'DisplayName', '\boldmath$\rho_{33} \ \mathrm{(Excited)}$');
grid on;
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Populations}$');
title('\boldmath$\mathrm{State \ Populations \ (Flat \ Signal)}$');
legend('show');
PlotUtils.styleAxes(gca);

% --- Figure 4: Zoom on Flat Excited State Population ---
figure('Color', 'w', 'Name', 'Diabatic Transition Zoom (Flat Signal)');
hold on;
plot(tgrid_Flat, rho_Flat(:,3), 'k', 'LineWidth', 4, 'DisplayName', '\boldmath$\rho_{33} \ \mathrm{(Excited)}$');
grid on;
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Population \ \rho_{33}}$');
title('\boldmath$\mathrm{Excited \ State \ Population \ (Flat \ Signal)}$');
xlim([245, 255]); % Zoom in precisely on the same region
legend('show');
PlotUtils.styleAxes(gca);

% =========================================================================
% 3. SIGNAL ANALYSIS: TIME DOMAIN, FFT, AND INSTANTANEOUS FREQUENCY
% =========================================================================

% Sampling parameters
f_sampling = 60 / (2 * pi);
dt = 1 / f_sampling;
T_period = compute_fundamental_period(angular_freqs_SO, f_sampling);
signals_duration = T_period * 100;
t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1)';

sampled_SO = SO_signal(t_axis);
sampled_Flat = Flat_signal(t_axis);

% --- Figure 5: Signals in Time Domain ---
figure('Color', 'w', 'Name', 'Signals in Time Domain');
hold on;
plot(t_axis, real(sampled_SO), 'r-', 'LineWidth', 2, 'DisplayName', '\textbf{SO}');
plot(t_axis, real(sampled_Flat), '-', 'Color', [0, 0.5, 0], 'LineWidth', 2, 'DisplayName', '\textbf{Flat}');
grid on;
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathrm{Amplitude \ [arb]}$');
title('\boldmath$\mathrm{Signals \ in \ Time \ Domain}$');
legend('show');
PlotUtils.styleAxes(gca);

% --- Figure 6: FFT of Signals ---
N_fft = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N_fft) * 2*pi;
fft_SO = fftshift(abs(fft(sampled_SO, N_fft))) / N_fft;
fft_Flat = fftshift(abs(fft(sampled_Flat, N_fft))) / N_fft;

figure('Color', 'w', 'Name', 'FFT of Signals');
hold on;
plot(freq_axis, fft_SO, 'r-', 'LineWidth', 2, 'DisplayName', '\textbf{SO}');
plot(freq_axis, fft_Flat, '-', 'Color', [0, 0.5, 0], 'LineWidth', 2, 'DisplayName', '\textbf{Flat}');
grid on;
xlabel('\boldmath$\mathrm{Angular \ frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathrm{Amplitude \ [arb]}$');
title('\boldmath$\mathrm{Signals \ in \ Frequency \ Domain \ (FFT)}$');
legend('show');
PlotUtils.styleAxes(gca);

% --- Figure 7: Instantaneous Frequency Analysis ---
t_inst = linspace(200, 300, 5000)';
dt_inst = t_inst(2) - t_inst(1);

y_so_inst   = arrayfun(SO_signal, t_inst);
y_flat_inst = arrayfun(Flat_signal, t_inst);

inst_freq_so = compute_instantaneous_frequency(y_so_inst, dt_inst);
inst_freq_flat = compute_instantaneous_frequency(y_flat_inst, dt_inst);

figure('Color', 'w', 'Name', 'Instantaneous Frequency Analysis');

x_limits = [245, 255];
y_limits = [-2, 4];

% Subplot 1: SO Instantaneous Frequency
subplot(2,1,1);
hold on;
plot(t_inst, real(y_so_inst), 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Wave');
plot(t_inst, inst_freq_so, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('SO: Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);

% Subplot 2: Flat Instantaneous Frequency
subplot(2,1,2);
hold on;
plot(t_inst, real(y_flat_inst), 'LineWidth', 2, 'Color', [0, 0.5, 0], 'DisplayName', 'Wave');
plot(t_inst, inst_freq_flat, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
title('Flat Spectrum: Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
PlotUtils.styleAxes(gca);
