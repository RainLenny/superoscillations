% TLS Dynamics Simulation: Optimization of 5-Harmonic Signal
clear; clc; close all;

%% 1. Parameters & Scaling
% --- ENERGY SCALING OPTION ---
Energy_Scale = 0.5; % Scale > 1 increases peak energy, < 1 decreases it.

% TLS parameters
T1 = 500;
T2 = 300;
Omega = 0.01;
rho30 = -1;

% Pulse parameters
T_pulse = 100;
t0 = 200;
A_SO = 7 * Energy_Scale;
A_ref1 = 25 * Energy_Scale;

% Frequencies
n = 1:5;
omega = 0.18 * n;

%% 2. Define Original Coefficients
% SO Coefficients from the paper
c_SO = [-0.156 + 0.3311i, ...
    -0.862 - 1.042i, ...
    2.341 - 0.601i, ...
    -0.502 + 2.634i, ...
    -1.820 - 1.322i];

% --- Find Target Peak Amplitude ---
t_search = linspace(0, 600, 10000);
SO_carrier_search = zeros(size(t_search));
for k = 1:length(t_search)
    SO_carrier_search(k) = sum(real(c_SO .* exp(1i * omega .* t_search(k))));
end
peak_SO_carrier = max(abs(SO_carrier_search));
target_peak_carrier = peak_SO_carrier; % We will constrain our optimizer to this

% --- Flat Spectrum Coefficients ---
peak_flat_carrier = 5; 
c_flat = (peak_SO_carrier / peak_flat_carrier) * ones(1, 5);

%% 3. Optimize Coefficients for Maximal Excitation
fprintf('Starting optimization. This may take a minute...\n');

% Initial guess: we start with the flat spectrum coefficients
x0 = [real(c_flat), imag(c_flat)];

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1500, ...
    'TolCon', 1e-4, ...
    'TolX', 1e-4);

% Run fmincon
% Objective: maximize max(rho3) -> minimize -max(rho3)
% Constraint: max(abs(carrier)) == target_peak_carrier
x_opt = fmincon(@(x) obj_fun(x, omega, T1, T2, Omega, rho30, t_search, T_pulse, t0, A_SO), ...
                x0, [], [], [], [], [], [], ...
                @(x) nonlcon(x, omega, t_search, target_peak_carrier), options);

% Reconstruct optimized complex coefficients
c_opt = x_opt(1:5) + 1i * x_opt(6:10);
fprintf('Optimization complete.\n\n');

%% 4. Define Signal Functions
% Original Signals
f_SO = @(t) A_SO * sum(real(c_SO .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);
f_ref1 = @(t) A_ref1 * cos(omega(5) * t) * exp(-(t-t0)^2 / T_pulse^2);
f_ref2 = @(t) A_SO * sum(real(c_flat .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);

% Optimized Signal
f_opt = @(t) A_SO * sum(real(c_opt .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);

%% 5. Generate Time-Domain Signals & Phase Dynamics
t_vec = linspace(0, 600, 5000);
sig_SO = zeros(size(t_vec)); sig_ref1 = zeros(size(t_vec));
sig_ref2 = zeros(size(t_vec)); sig_opt = zeros(size(t_vec));
d_angle_SO = zeros(size(t_vec)); d_angle_ref2 = zeros(size(t_vec));
d_angle_opt = zeros(size(t_vec));

for k = 1:length(t_vec)
    % Amplitudes
    sig_SO(k) = f_SO(t_vec(k));
    sig_ref1(k) = f_ref1(t_vec(k));
    sig_ref2(k) = f_ref2(t_vec(k));
    sig_opt(k) = f_opt(t_vec(k));
    
    % Instantaneous Frequencies
    Z_SO = sum(c_SO .* exp(1i * omega .* t_vec(k)));
    dZ_SO = sum(1i * omega .* c_SO .* exp(1i * omega .* t_vec(k)));
    d_angle_SO(k) = imag(dZ_SO / Z_SO);
    
    Z_ref2 = sum(c_flat .* exp(1i * omega .* t_vec(k)));
    dZ_ref2 = sum(1i * omega .* c_flat .* exp(1i * omega .* t_vec(k)));
    d_angle_ref2(k) = imag(dZ_ref2 / Z_ref2);
    
    Z_opt = sum(c_opt .* exp(1i * omega .* t_vec(k)));
    dZ_opt = sum(1i * omega .* c_opt .* exp(1i * omega .* t_vec(k)));
    d_angle_opt(k) = imag(dZ_opt / Z_opt);
end

%% 6. Calculate Frequency-Domain Signals (Spectral Density)
nu_vec = linspace(0, 1.2, 1000);
spec_SO = zeros(size(nu_vec)); spec_ref1 = zeros(size(nu_vec));
spec_ref2 = zeros(size(nu_vec)); spec_opt = zeros(size(nu_vec));

for k = 1:length(nu_vec)
    exp_term = exp(1i * nu_vec(k) * t_vec);
    spec_SO(k) = abs(trapz(t_vec, sig_SO .* exp_term));
    spec_ref1(k) = abs(trapz(t_vec, sig_ref1 .* exp_term));
    spec_ref2(k) = abs(trapz(t_vec, sig_ref2 .* exp_term));
    spec_opt(k) = abs(trapz(t_vec, sig_opt .* exp_term));
end

%% 7. Solve Density Matrix ODEs (Eq 5)
rho_init = [0; 0; -1];
tspan = [0 600];

ode_system = @(t, rho, func) [ ...
    rho(2) - rho(1)/T2; ...
    -rho(1) - rho(2)/T2 + 2*Omega*func(t)*rho(3); ...
    -2*Omega*func(t)*rho(2) - (rho(3)-rho30)/T1 ];

[t_SO, rho_SO] = ode45(@(t, rho) ode_system(t, rho, f_SO), tspan, rho_init);
[t_ref1, rho_ref1] = ode45(@(t, rho) ode_system(t, rho, f_ref1), tspan, rho_init);
[t_ref2, rho_ref2] = ode45(@(t, rho) ode_system(t, rho, f_ref2), tspan, rho_init);
[t_opt, rho_opt] = ode45(@(t, rho) ode_system(t, rho, f_opt), tspan, rho_init);

%% 8. Plotting
% --- Figure 1: Signals in Time ---
figure('Name', 'Signals in Time', 'Position', [100, 100, 800, 800]);

subplot(3, 1, 1);
yyaxis left; plot(t_vec, Omega * sig_SO, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); ylabel('\Omega f_{SO}(t)');
yyaxis right; plot(t_vec, d_angle_SO, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5); ylabel('d(angle)/dt');
title('Superoscillating (SO) Signal'); xlim([0 600]); grid on;

subplot(3, 1, 2);
yyaxis left; plot(t_vec, Omega * sig_ref2, '-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5); ylabel('\Omega f_{ref2}(t)');
yyaxis right; plot(t_vec, d_angle_ref2, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5); ylabel('d(angle)/dt');
title('Ref 2: Flat Spectrum'); xlim([0 600]); grid on;

subplot(3, 1, 3);
yyaxis left; plot(t_vec, Omega * sig_opt, '-', 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1.5); ylabel('\Omega f_{opt}(t)');
yyaxis right; plot(t_vec, d_angle_opt, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5); ylabel('d(angle)/dt');
title('Optimized Spectrum (Maximized Excitation)'); xlabel('Time'); xlim([0 600]); grid on;

% --- Figure 2: Signals in Frequency ---
figure('Name', 'Signals in Frequency');
plot(nu_vec, spec_SO, 'b', 'LineWidth', 1); hold on;
plot(nu_vec, spec_ref1, 'r--', 'LineWidth', 1);
plot(nu_vec, spec_ref2, 'g:', 'LineWidth', 1.5);
plot(nu_vec, spec_opt, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 2);
xline(1.0, 'k-', 'Absorption Band (\omega_0=1)');
title('Signals in Frequency (Spectral Density)');
xlabel('Frequency'); ylabel('Amplitude (a.u.)');
xlim([0 1.2]); grid on;
legend('SO', 'Single Harmonic', 'Flat', 'Optimized', 'Location', 'best');

% --- Figure 3: Excitation in Time ---
figure('Name', 'Excitation in Time');
plot(t_SO, rho_SO(:,3), 'b', 'LineWidth', 1.5); hold on;
plot(t_ref1, rho_ref1(:,3), 'r--', 'LineWidth', 1.5);
plot(t_ref2, rho_ref2(:,3), 'g:', 'LineWidth', 1.5);
plot(t_opt, rho_opt(:,3), 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 2);
yline(0, 'k--');
title('Excitation in Time (\rho_3 = \rho_{22} - \rho_{11})');
xlabel('Time'); ylabel('Population Inversion (\rho_3)');
ylim([-1 1]); xlim([0 600]); grid on;
legend('SO', 'Single Harmonic', 'Flat', 'Optimized', 'Location', 'best');


%% LOCAL FUNCTIONS FOR OPTIMIZATION
function neg_rho3_max = obj_fun(x, omega, T1, T2, Omega, rho30, t_search, T_pulse, t0, A_SO)
    % Reconstruct complex coefficients
    c = x(1:5) + 1i * x(6:10);
    
    % Reconstruct signal
    func = @(t) A_SO * sum(real(c .* exp(1i * omega .* t))) .* exp(-(t-t0).^2 / T_pulse^2);
    
    % Solve ODE
    ode_sys = @(t, rho) [ ...
        rho(2) - rho(1)/T2; ...
        -rho(1) - rho(2)/T2 + 2*Omega*func(t)*rho(3); ...
        -2*Omega*func(t)*rho(2) - (rho(3)-rho30)/T1 ];
    
    % Fast evaluate to find max excitation
    [~, rho] = ode45(ode_sys, [0 600], [0; 0; -1]);
    
    % We want to maximize the maximum population inversion
    neg_rho3_max = -max(rho(:,3)); 
end

function [c_ineq, ceq] = nonlcon(x, omega, t_search, target_peak)
    % Reconstruct complex coefficients
    c = x(1:5) + 1i * x(6:10);
    
    % Calculate carrier wave over search space
    carrier = zeros(size(t_search));
    for k = 1:length(t_search)
        carrier(k) = sum(real(c .* exp(1i * omega .* t_search(k))));
    end
    
    % Enforce that peak carrier amplitude matches the target exactly
    ceq = max(abs(carrier)) - target_peak; 
    c_ineq = [];
end