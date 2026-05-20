% TLS Dynamics Simulation: Verifying Excitation Mechanism via Peak Isolation
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

%% 2. Define Signal Functions & Thresholding
% Original Signals
f_SO = @(t) A_SO * sum(real(c_SO .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);
f_ref1 = @(t) A_ref1 * cos(omega(5) * t) * exp(-(t-t0)^2 / T_pulse^2);
f_ref2 = @(t) A_SO * sum(real(c_flat .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);

% --- Mechanism Verification: Peak Isolation (Thresholding) ---
% Set the threshold just above the amplitude of the superoscillating ripples
threshold = 40; 

% Hard threshold: Set to 0 if the absolute amplitude is below the threshold.
% This zeros out the superoscillatory region but leaves the high peaks intact.
f_SO_peaks = @(t) f_SO(t) * (abs(f_SO(t)) >= threshold);
f_ref2_peaks = @(t) f_ref2(t) * (abs(f_ref2(t)) >= threshold);

%% 3. Generate Time-Domain Signals for Plotting
t_vec = linspace(0, 600, 5000);
sig_SO = zeros(size(t_vec));
sig_ref2 = zeros(size(t_vec));
sig_SO_peaks_vec = zeros(size(t_vec));
sig_ref2_peaks_vec = zeros(size(t_vec));

for k = 1:length(t_vec)
    sig_SO(k) = f_SO(t_vec(k));
    sig_ref2(k) = f_ref2(t_vec(k));
    sig_SO_peaks_vec(k) = f_SO_peaks(t_vec(k));
    sig_ref2_peaks_vec(k) = f_ref2_peaks(t_vec(k));
end

%% 4. Solve Density Matrix ODEs (Eq 5)
% Initial conditions: rho1=0, rho2=0, rho3=-1
rho_init = [0; 0; -1];
tspan = [0 600];

% ODE system function for the Two-Level System
ode_system = @(t, rho, func) [ ...
    rho(2) - rho(1)/T2; ...
    -rho(1) - rho(2)/T2 + 2*Omega*func(t)*rho(3); ...
    -2*Omega*func(t)*rho(2) - (rho(3)-rho30)/T1 ];

% Solve for original signals
[t_SO, rho_SO] = ode45(@(t, rho) ode_system(t, rho, f_SO), tspan, rho_init);
[t_ref2, rho_ref2] = ode45(@(t, rho) ode_system(t, rho, f_ref2), tspan, rho_init);

% Solve for peak-only signals
[t_SO_peaks, rho_SO_peaks] = ode45(@(t, rho) ode_system(t, rho, f_SO_peaks), tspan, rho_init);
[t_ref2_peaks, rho_ref2_peaks] = ode45(@(t, rho) ode_system(t, rho, f_ref2_peaks), tspan, rho_init);

%% 5. Plotting: Verification of Excitation Mechanism
figure('Name', 'Mechanism Verification: Peak-Only Signals', 'Position', [100, 100, 900, 700]);

% Subplot 1: Compare Original vs. Peak-Only Signals in Time
ax1 = subplot(2, 1, 1); 
plot(t_vec, sig_SO, 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 2); hold on;
plot(t_vec, sig_SO_peaks_vec, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
plot(t_vec, sig_ref2, 'Color', [0.4660 0.6740 0.1880 0.3], 'LineWidth', 2);
plot(t_vec, sig_ref2_peaks_vec, ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);

% Plot threshold boundaries
yline(threshold, 'k--', 'Threshold Limit', 'LabelHorizontalAlignment', 'left');
yline(-threshold, 'k--');

title('Time Domain: Original vs. Peak-Isolated Signals');
xlabel('Time'); 
ylabel('Amplitude');
xlim([0 600]); 
grid on;
legend('SO Original', 'SO Peaks Only', 'Ref2 Original', 'Ref2 Peaks Only', 'Location', 'best');

% Subplot 2: Resulting Population Inversion
ax2 = subplot(2, 1, 2); 
plot(t_SO, rho_SO(:,3), 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 2); hold on;
plot(t_SO_peaks, rho_SO_peaks(:,3), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
plot(t_ref2, rho_ref2(:,3), 'Color', [0.4660 0.6740 0.1880 0.3], 'LineWidth', 2);
plot(t_ref2_peaks, rho_ref2_peaks(:,3), ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);

% Zero inversion threshold
yline(0, 'k-', 'Zero Inversion', 'LabelHorizontalAlignment', 'left');

title('Excitation Verification (\rho_3 = \rho_{22} - \rho_{11})');
xlabel('Time'); 
ylabel('Population Inversion (\rho_3)');
ylim([-1 1.1]); 
xlim([0 600]); 
grid on;
legend('SO Full Excitation', 'SO Peaks Excitation', 'Ref2 Full Excitation', 'Ref2 Peaks Excitation', 'Location', 'best');

% Link the X-axes together
linkaxes([ax1, ax2], 'x');