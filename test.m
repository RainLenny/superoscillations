% TLS Dynamics Simulation: Superoscillating Field vs. References
clear; clc; close all;

%% 1. Parameters
% TLS parameters
T1 = 500;       
T2 = 300;       
Omega = 0.01;   
rho30 = -1;     

% Pulse parameters
T_pulse = 100;
t0 = 200;
A_SO = 7;
A_ref1 = 25;

% Frequencies
n = 1:5;
omega = 0.18 * n; 

% SO Coefficients from the paper
c_SO = [-0.156 + 0.3311i, ...
        -0.862 - 1.042i, ...
         2.341 - 0.601i, ...
        -0.502 + 2.634i, ...
        -1.820 - 1.322i];

% Flat Spectrum Coefficients (Same total energy, equally distributed)
c_rms = sqrt(mean(abs(c_SO).^2));
c_flat = c_rms * ones(1, 5);

%% 2. Define Signal Functions
% Superoscillating signal
f_SO = @(t) A_SO * sum(real(c_SO .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);

% Reference 1: Fastest harmonic only (omega_5 = 0.90)
f_ref1 = @(t) A_ref1 * cos(omega(5) * t) * exp(-(t-t0)^2 / T_pulse^2);

% Reference 2: Flat spectrum (equal energy per frequency)
f_ref2 = @(t) A_SO * sum(real(c_flat .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);

%% 3. Generate Time-Domain Signals
t_vec = linspace(0, 600, 5000);
sig_SO = zeros(size(t_vec));
sig_ref1 = zeros(size(t_vec));
sig_ref2 = zeros(size(t_vec));

for k = 1:length(t_vec)
    sig_SO(k) = f_SO(t_vec(k));
    sig_ref1(k) = f_ref1(t_vec(k));
    sig_ref2(k) = f_ref2(t_vec(k));
end

%% 4. Calculate Frequency-Domain Signals (Spectral Density)
nu_vec = linspace(0, 1.2, 1000);
spec_SO = zeros(size(nu_vec));
spec_ref1 = zeros(size(nu_vec));
spec_ref2 = zeros(size(nu_vec));

% Numerical integration for Fourier transform: abs(int(f(t)*exp(i*nu*t) dt))
for k = 1:length(nu_vec)
    exp_term = exp(1i * nu_vec(k) * t_vec);
    spec_SO(k) = abs(trapz(t_vec, sig_SO .* exp_term));
    spec_ref1(k) = abs(trapz(t_vec, sig_ref1 .* exp_term));
    spec_ref2(k) = abs(trapz(t_vec, sig_ref2 .* exp_term));
end

%% 5. Solve Density Matrix ODEs (Eq 5)
% Initial conditions: rho1=0, rho2=0, rho3=-1
rho_init = [0; 0; -1];
tspan = [0 600];

% ODE system function
ode_system = @(t, rho, func) [ ...
    rho(2) - rho(1)/T2; ...
    -rho(1) - rho(2)/T2 + 2*Omega*func(t)*rho(3); ...
    -2*Omega*func(t)*rho(2) - (rho(3)-rho30)/T1 ];

% Solve for each signal
[t_SO, rho_SO] = ode45(@(t, rho) ode_system(t, rho, f_SO), tspan, rho_init);
[t_ref1, rho_ref1] = ode45(@(t, rho) ode_system(t, rho, f_ref1), tspan, rho_init);
[t_ref2, rho_ref2] = ode45(@(t, rho) ode_system(t, rho, f_ref2), tspan, rho_init);

%% 6. Plotting
figure('Position', [100, 100, 1200, 400]);

% --- Subplot 1: Signals in Time ---
subplot(1, 3, 1);
plot(t_vec, Omega * sig_SO, 'b', 'LineWidth', 1); hold on;
plot(t_vec, Omega * sig_ref1, 'r--', 'LineWidth', 1);
plot(t_vec, Omega * sig_ref2, 'g:', 'LineWidth', 1.5);
title('Signals in Time (Rabi Frequency \Omega f(t))');
xlabel('Time');
ylabel('\Omega f(t)');
legend('Superoscillating (SO)', 'Ref 1: 0.9\omega_0', 'Ref 2: Flat Spectrum', 'Location', 'best');
xlim([0 600]);
grid on;

% --- Subplot 2: Signals in Frequency ---
subplot(1, 3, 2);
plot(nu_vec, spec_SO, 'b', 'LineWidth', 1); hold on;
plot(nu_vec, spec_ref1, 'r--', 'LineWidth', 1);
plot(nu_vec, spec_ref2, 'g:', 'LineWidth', 1.5);
xline(1.0, 'k-', 'Absorption Band (\omega_0=1)'); % Marker for TLS frequency
title('Signals in Frequency (Spectral Density)');
xlabel('Frequency');
ylabel('Amplitude (a.u.)');
xlim([0 1.2]);
grid on;

% --- Subplot 3: Excitation in Time ---
subplot(1, 3, 3);
plot(t_SO, rho_SO(:,3), 'b', 'LineWidth', 1.5); hold on;
plot(t_ref1, rho_ref1(:,3), 'r--', 'LineWidth', 1.5);
plot(t_ref2, rho_ref2(:,3), 'g:', 'LineWidth', 1.5);
yline(0, 'k--'); % Marker for positive inversion threshold
title('Excitation in Time (\rho_3 = \rho_{22} - \rho_{11})');
xlabel('Time');
ylabel('Population Inversion (\rho_3)');
ylim([-1 1]);
xlim([0 600]);
grid on;