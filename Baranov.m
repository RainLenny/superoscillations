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

%% 3.5 Calculate Instantaneous Frequency (d_angle/dt) Analytically
% The analytical derivative of the phase is given by Im(Z'(t) / Z(t))
% where Z(t) is the complex carrier signal.

d_angle_SO = zeros(size(t_vec));
d_angle_ref1 = zeros(size(t_vec));
d_angle_ref2 = zeros(size(t_vec));

for k = 1:length(t_vec)
    % Superoscillating (SO)
    Z_SO_val = sum(c_SO .* exp(1i * omega .* t_vec(k)));
    dZ_SO_val = sum(1i * omega .* c_SO .* exp(1i * omega .* t_vec(k)));
    d_angle_SO(k) = imag(dZ_SO_val / Z_SO_val);

    % Ref 1 (Single Harmonic)
    Z_ref1_val = exp(1i * omega(5) * t_vec(k));
    dZ_ref1_val = 1i * omega(5) * exp(1i * omega(5) * t_vec(k));
    d_angle_ref1(k) = imag(dZ_ref1_val / Z_ref1_val);

    % Ref 2 (Flat Spectrum)
    Z_ref2_val = sum(c_flat .* exp(1i * omega .* t_vec(k)));
    dZ_ref2_val = sum(1i * omega .* c_flat .* exp(1i * omega .* t_vec(k)));
    d_angle_ref2(k) = imag(dZ_ref2_val / Z_ref2_val);
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
% --- Figure 1: Signals in Time ---
figure('Name', 'Signals in Time', 'Position', [100, 100, 800, 800]);

% Subplot 1: Superoscillating (SO)
subplot(3, 1, 1);
yyaxis left
plot(t_vec, Omega * sig_SO, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
ylabel('\Omega f_{SO}(t)');
yyaxis right
plot(t_vec, d_angle_SO, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('d(angle)/dt');
title('Superoscillating (SO) Signal');
xlim([0 600]); grid on;
legend('Wave', 'd\_angle/dt', 'Location', 'best');

% Subplot 2: Reference 1
subplot(3, 1, 3);
yyaxis left
plot(t_vec, Omega * sig_ref1, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('\Omega f_{ref1}(t)');
yyaxis right
plot(t_vec, d_angle_ref1, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5);
ylabel('d(angle)/dt');
title('Ref 1: Fastest Harmonic (0.9\omega_0)');
xlim([0 600]);
ylim([0 1.5]); % Manually set so the flat line isn't zoomed in infinitely
grid on;
legend('Wave', 'd\_angle/dt', 'Location', 'best');

% Subplot 3: Reference 2
subplot(3, 1, 2);
yyaxis left
plot(t_vec, Omega * sig_ref2, '-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
ylabel('\Omega f_{ref2}(t)');
yyaxis right
plot(t_vec, d_angle_ref2, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('d(angle)/dt');
title('Ref 2: Flat Spectrum');
xlabel('Time'); xlim([0 600]); grid on;
legend('Wave', 'd\_angle/dt', 'Location', 'best');

% --- Figure 2: Signals in Frequency ---
figure('Name', 'Signals in Frequency');
plot(nu_vec, spec_SO, 'b', 'LineWidth', 1); hold on;
plot(nu_vec, spec_ref1, 'r--', 'LineWidth', 1);
plot(nu_vec, spec_ref2, 'g:', 'LineWidth', 1.5);
xline(1.0, 'k-', 'Absorption Band (\omega_0=1)'); % Marker for TLS frequency
title('Signals in Frequency (Spectral Density)');
xlabel('Frequency');
ylabel('Amplitude (a.u.)');
xlim([0 1.2]);
grid on;

% --- Figure 2: Excitation in Time ---
figure('Name', 'Excitation in Time');
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

% --- Figure 4: Signals in Time (No Envelope) ---
% First, compute the continuous signals without the exp(-(t-t0)^2 / T_pulse^2) term
sig_SO_cont = zeros(size(t_vec));
sig_ref1_cont = zeros(size(t_vec));
sig_ref2_cont = zeros(size(t_vec));

for k = 1:length(t_vec)
    sig_SO_cont(k) = A_SO * sum(real(c_SO .* exp(1i * omega .* t_vec(k))));
    sig_ref1_cont(k) = A_ref1 * cos(omega(5) * t_vec(k));
    sig_ref2_cont(k) = A_SO * sum(real(c_flat .* exp(1i * omega .* t_vec(k))));
end

% Now plot them just like Figure 1
figure('Name', 'Continuous Signals in Time (No Envelope)', 'Position', [200, 200, 800, 800]);

% Subplot 1: Superoscillating (SO) Without Envelope
subplot(3, 1, 1);
yyaxis left
plot(t_vec, Omega * sig_SO_cont, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
ylabel('\Omega f_{SO, cont}(t)');
yyaxis right
plot(t_vec, d_angle_SO, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('d(angle)/dt');
title('Superoscillating (SO) Signal (No Envelope)');
xlim([0 600]); grid on;
legend('Continuous Wave', 'd\_angle/dt', 'Location', 'best');

% Subplot 2: Reference 1 Without Envelope
subplot(3, 1, 3);
yyaxis left
plot(t_vec, Omega * sig_ref1_cont, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('\Omega f_{ref1, cont}(t)');
yyaxis right
plot(t_vec, d_angle_ref1, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5);
ylabel('d(angle)/dt');
title('Ref 1: Fastest Harmonic (No Envelope)');
xlim([0 600]);
ylim([0 1.5]); % Constrain flat d_angle/dt to prevent infinite zooming
grid on;
legend('Continuous Wave', 'd\_angle/dt', 'Location', 'best');

% Subplot 3: Reference 2 Without Envelope
subplot(3, 1, 2);
yyaxis left
plot(t_vec, Omega * sig_ref2_cont, '-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
ylabel('\Omega f_{ref2, cont}(t)');
yyaxis right
plot(t_vec, d_angle_ref2, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('d(angle)/dt');
title('Ref 2: Flat Spectrum (No Envelope)');
xlabel('Time'); xlim([0 600]); grid on;
legend('Continuous Wave', 'd\_angle/dt', 'Location', 'best');

% %% other plot 
% figure
% hold on
% 
% plot(t_vec, Omega * sig_SO, 'b', 'LineWidth', 1);
% plot(t_vec, Omega * sig_ref2, 'g', 'LineWidth', 1.5);
% 
% figure
% hold on
% 
% plot(nu_vec, spec_SO, 'b', 'LineWidth', 1);
% plot(nu_vec, spec_ref2, 'g', 'LineWidth', 1.5);
