% TLS Dynamics Simulation: Superoscillating Field vs. References
clear; clc; close all;

%% 1. Parameters & Global Scaling
% --- Global Amplitude Scale (Peak Rabi Frequency Control) ---
global_scaling = 40;  % Single parameter to scale all signal peaks simultaneously

% TLS parameters
T1 = 500;
T2 = 300;
Omega = 0.01;
rho30 = -1;

% Pulse parameters
T_pulse = 100;
t0 = 200;

% Frequencies & Coefficients (Defined as column vectors for vectorization)
omega = 0.18 * (1:5)'; 
c_SO = [-0.156 + 0.3311i; ...
        -0.862 - 1.042i; ...
         2.341 - 0.601i; ...
        -0.502 + 2.634i; ...
        -1.820 - 1.322i];

c_flat = ones(5, 1);

%% 2. Generate Time-Domain Signals & Instantaneous Frequencies
t_vec = linspace(0, 600, 5000);
envelope = exp(-(t_vec - t0).^2 / T_pulse^2);

% --- Vectorized Complex Carriers & Derivatives ---
Z_SO = sum(c_SO .* exp(1i * omega * t_vec), 1);
dZ_SO = sum((1i * omega .* c_SO) .* exp(1i * omega * t_vec), 1);

Z_ref2 = sum(c_flat .* exp(1i * omega * t_vec), 1);
dZ_ref2 = sum((1i * omega .* c_flat) .* exp(1i * omega * t_vec), 1);

% --- Raw (Unscaled) Time Signals ---
sig_SO_raw   = real(Z_SO) .* envelope;
sig_ref1_raw = cos(omega(5) * t_vec) .* envelope; % Cosine reference (Fastest harmonic)
sig_ref2_raw = real(Z_ref2) .* envelope;

% --- Analytical Instantaneous Frequencies (d_angle/dt = Im(Z'/Z)) ---
d_angle_SO   = imag(dZ_SO ./ Z_SO);
d_angle_ref1 = omega(5) * ones(size(t_vec)); % Pure harmonic is constant
d_angle_ref2 = imag(dZ_ref2 ./ Z_ref2);

% --- Exact Peak Normalization ---
% Finds the absolute peak of each raw signal and normalizes it to A_global
N_SO   = max(abs(sig_SO_raw));
N_ref1 = max(abs(sig_ref1_raw));
N_ref2 = max(abs(sig_ref2_raw));

sig_SO   = global_scaling * (sig_SO_raw / N_SO);
sig_ref1 = global_scaling * (sig_ref1_raw / N_ref1);
sig_ref2 = global_scaling * (sig_ref2_raw / N_ref2);

%% 3. Calculate Frequency-Domain Signals (Spectral Density)
nu_vec = linspace(0, 1.2, 1000);
spec_SO = zeros(size(nu_vec));
spec_ref1 = zeros(size(nu_vec));
spec_ref2 = zeros(size(nu_vec));

% Numerical Fourier Transform
for k = 1:length(nu_vec)
    exp_term = exp(1i * nu_vec(k) * t_vec);
    spec_SO(k)   = abs(trapz(t_vec, sig_SO .* exp_term));
    spec_ref1(k) = abs(trapz(t_vec, sig_ref1 .* exp_term));
    spec_ref2(k) = abs(trapz(t_vec, sig_ref2 .* exp_term));
end

%% 4. Solve Density Matrix ODEs
rho_init = [0; 0; -1];
tspan = [0 600];

% Inline scalar helper functions for the ODE solver (retains exact normalization)
envelope_t    = @(t) exp(-(t-t0)^2 / T_pulse^2);
f_SO_scalar   = @(t) global_scaling * real(sum(c_SO .* exp(1i * omega * t))) * envelope_t(t) / N_SO;
f_ref1_scalar = @(t) global_scaling * cos(omega(5) * t) * envelope_t(t) / N_ref1;
f_ref2_scalar = @(t) global_scaling * real(sum(c_flat .* exp(1i * omega * t))) * envelope_t(t) / N_ref2;

ode_system = @(t, rho, func) [ ...
    rho(2) - rho(1)/T2; ...
   -rho(1) - rho(2)/T2 + 2*Omega*func(t)*rho(3); ...
   -2*Omega*func(t)*rho(2) - (rho(3)-rho30)/T1 ];

[t_SO, rho_SO]     = ode45(@(t, rho) ode_system(t, rho, f_SO_scalar), tspan, rho_init);
[t_ref1, rho_ref1] = ode45(@(t, rho) ode_system(t, rho, f_ref1_scalar), tspan, rho_init);
[t_ref2, rho_ref2] = ode45(@(t, rho) ode_system(t, rho, f_ref2_scalar), tspan, rho_init);

%% 5. Plotting
% --- Figure 1: Signals in Time (3 Subplots) ---
figure('Name', 'Signals in Time', 'Position', [100, 100, 800, 600]);
freq_color = [0.4940 0.1840 0.5560]; % Consistent purple color for frequency axis

% Subplot 1: Superoscillating (SO)
subplot(3, 1, 1);
yyaxis left
plot(t_vec, Omega * sig_SO, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
ylabel('\Omega f_{SO}(t)'); ylim([-global_scaling*Omega*1.1, global_scaling*Omega*1.1]);
yyaxis right
plot(t_vec, d_angle_SO, '-', 'Color', freq_color, 'LineWidth', 1.5);
ylabel('d(angle)/dt'); ylim([0 3]);
title('Superoscillating (SO) Signal'); xlim([0 600]); grid on;

% Subplot 2: Reference 1 (Single Harmonic Cosine)
subplot(3, 1, 2);
yyaxis left
plot(t_vec, Omega * sig_ref1, '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.5);
ylabel('\Omega f_{ref1}(t)'); ylim([-global_scaling*Omega*1.1, global_scaling*Omega*1.1]);
yyaxis right
plot(t_vec, d_angle_ref1, '-', 'Color', freq_color, 'LineWidth', 1.5);
ylabel('d(angle)/dt'); ylim([0 3]);
title('Ref 1: Single Harmonic (\omega_5)'); xlim([0 600]); grid on;

% Subplot 3: Reference 2 (Flat Spectrum)
subplot(3, 1, 3);
yyaxis left
plot(t_vec, Omega * sig_ref2, '-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
ylabel('\Omega f_{ref2}(t)'); ylim([-global_scaling*Omega*1.1, global_scaling*Omega*1.1]);
yyaxis right
plot(t_vec, d_angle_ref2, '-', 'Color', freq_color, 'LineWidth', 1.5);
ylabel('d(angle)/dt'); ylim([0 3]);
title('Ref 2: Flat Spectrum'); xlabel('Time'); xlim([0 600]); grid on;

% --- Figure 2: Signals in Frequency ---
figure('Name', 'Signals in Frequency');
plot(nu_vec, spec_SO, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
plot(nu_vec, spec_ref1, 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 1.5);
plot(nu_vec, spec_ref2, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':', 'LineWidth', 2);
xline(1.0, 'k-', 'Absorption Band (\omega_0=1)', 'LabelHorizontalAlignment', 'center'); 
title('Signals in Frequency (Spectral Density)');
xlabel('Frequency'); ylabel('Amplitude (a.u.)');
xlim([0 1.2]); legend('SO', 'Ref 1 (Harmonic)', 'Ref 2 (Flat)', 'Location', 'best'); grid on;

% --- Figure 3: Excitation in Time ---
figure('Name', 'Excitation in Time');
plot(t_SO, rho_SO(:,3), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
plot(t_ref1, rho_ref1(:,3), 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 1.5);
plot(t_ref2, rho_ref2(:,3), 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':', 'LineWidth', 2);
yline(0, 'k--'); 
title('Excitation in Time (\rho_3 = \rho_{22} - \rho_{11})');
xlabel('Time'); ylabel('Population Inversion (\rho_3)');
ylim([-1 1]); xlim([0 600]);
legend('SO', 'Ref 1 (Harmonic)', 'Ref 2 (Flat)', 'Location', 'best'); grid on;