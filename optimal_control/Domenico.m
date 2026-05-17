% Optimal Control of Superoscillating Excitation in a Two-Level System
% Based on Baranov et al. (2014) and D'Alessandro & Dahleh (2001) https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=928587&casa_token=fLaZXv6dlMoAAAAA:Qp50qS5kg-RNG-wBi_YGTgMXi_8AozujNiHMtM6mpCYOOAey0nxmpoNFYGSAxPZkrDCzdApz&tag=1

clear; clc; close all;

%% 1. System and Field Parameters
w0 = 1;                 % TLS resonance frequency
T1 = 500;               % Spontaneous decay time
T2 = 300;               % Dephasing time
Omega = 0.01;           % Rabi constant
T_pulse = 100;          % Pulse duration
t0 = 200;               % Pulse center
t_span = [0 600];       % Simulation time span

% The 5 non-resonant frequencies (all < w0 = 1)
w = 0.18 * (1:5);       

% Original Baranov Coefficients (Eq. 6)
C_bar = [-0.156 + 0.3311i, ...
         -0.862 - 1.042i, ...
          2.341 - 0.601i, ...
         -0.502 + 2.634i, ...
         -1.820 - 1.322i];
A_bar = 7; % Base amplitude factor

%% 2. Original Baranov Signal Definition & Peak Amplitude
% We need the peak amplitude of the Baranov signal to normalize our optimal signal.
t_fine = linspace(t_span(1), t_span(2), 3000);
s_bar_fine = zeros(size(t_fine));
for n = 1:5
    s_bar_fine = s_bar_fine + real(C_bar(n) * exp(1i * w(n) * t_fine));
end
f_bar_fine = A_bar * s_bar_fine .* exp(-((t_fine - t0)/T_pulse).^2);
max_f_bar = max(abs(f_bar_fine)); % Peak Electric Field (constrains the optimal control)

%% 3. Optimal Control - Parameter Optimization
% We optimize the real and imaginary parts of the 5 frequency coefficients.
% Initial guess is the Baranov signal to ensure we find a strictly better local/global minimum.
initial_guess = [real(C_bar), imag(C_bar)];

disp('Running Optimization for the Control Field... (This may take a minute)');
options = optimset('MaxFunEvals', 3000, 'MaxIter', 1500, 'Display', 'iter');

% Objective function minimizes the ground state / maximizes excited state at t_end
optimal_vars = fminsearch(@(vars) ...
    cost_function(vars, w, t0, T_pulse, max_f_bar, Omega, T1, T2, t_span), ...
    initial_guess, options);

% Reconstruct optimal coefficients and normalize them to the constraint
C_opt_raw = optimal_vars(1:5) + 1i * optimal_vars(6:10);
s_opt_fine = zeros(size(t_fine));
for n = 1:5
    s_opt_fine = s_opt_fine + real(C_opt_raw(n) * exp(1i * w(n) * t_fine));
end
f_opt_raw = s_opt_fine .* exp(-((t_fine - t0)/T_pulse).^2);
max_f_opt = max(abs(f_opt_raw));

% Normalized Optimal Coefficients
C_opt = C_opt_raw * (max_f_bar / max_f_opt); 

%% 4. Solve Optical Bloch Equations for Both Signals
disp('Simulating Optical Bloch Equations...');

% Initial state: rho1 = 0, rho2 = 0, rho3 = -1 (Ground state)
y0 = [0; 0; -1];

% Solve Baranov
[t_bar, rho_bar] = ode45(@(t, y) obe_solver(t, y, C_bar*A_bar, w, t0, T_pulse, Omega, T1, T2), t_span, y0);

% Solve Optimal
[t_opt, rho_opt] = ode45(@(t, y) obe_solver(t, y, C_opt, w, t0, T_pulse, Omega, T1, T2), t_span, y0);

% Excitation probabilities: rho22 = (rho3 + 1)/2
P_exc_bar = (rho_bar(:,3) + 1) / 2;
P_exc_opt = (rho_opt(:,3) + 1) / 2;

%% 5. Generate Time & Frequency Domain Signals for Plotting
f_bar_plot = calc_field(t_fine, C_bar*A_bar, w, t0, T_pulse);
f_opt_plot = calc_field(t_fine, C_opt, w, t0, T_pulse);

% Compute FFT for frequency domain
dt = t_fine(2) - t_fine(1);
Fs = 1/dt;
N = length(t_fine);
freqs = linspace(-Fs/2, Fs/2, N) * 2 * pi;

F_bar = fftshift(abs(fft(f_bar_plot))) / N;
F_opt = fftshift(abs(fft(f_opt_plot))) / N;

%% 6. Plotting Results
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Time Domain Driving Signals
subplot(2, 2, 1);
plot(t_fine, Omega * f_bar_plot, 'b', 'LineWidth', 1.5); hold on;
plot(t_fine, Omega * f_opt_plot, 'r--', 'LineWidth', 1.5);
xlabel('Time (\omega_0^{-1})'); ylabel('\Omega f(t)');
title('Driving Electric Fields (Rabi Frequency)');
legend('Baranov Original', 'Optimal Control');
grid on; xlim([0 600]);

% Plot 2: Frequency Spectrum
subplot(2, 2, 2);
plot(freqs, F_bar, 'b', 'LineWidth', 1.5); hold on;
plot(freqs, F_opt, 'r--', 'LineWidth', 1.5);
xline(w0, 'k:', 'LineWidth', 2, 'Label', '\omega_0 (Resonance)');
xlabel('Frequency (\omega)'); ylabel('Spectral Amplitude');
title('Field Spectrum (All frequencies < \omega_0)');
legend('Baranov Original', 'Optimal Control', 'Location', 'NorthWest');
xlim([0 1.2]); grid on;

% Plot 3: Excitation Probability (rho22)
subplot(2, 2, [3, 4]);
plot(t_bar, P_exc_bar, 'b', 'LineWidth', 2); hold on;
plot(t_opt, P_exc_opt, 'r', 'LineWidth', 2);
xlabel('Time (\omega_0^{-1})'); ylabel('Excitation Probability \rho_{22}');
title('TLS Excitation Dynamics');
legend('Baranov Original', sprintf('Optimal Control (Final Excitation: %.3f)', P_exc_opt(end)));
grid on; xlim([0 600]);

%% Helper Functions

% Calculates the objective cost function (Maximizing rho3 at the end of the pulse)
function cost = cost_function(vars, w, t0, T_pulse, max_f_bar, Omega, T1, T2, t_span)
    C_raw = vars(1:5) + 1i * vars(6:10);
    
    % Ensure peak constraint normalization
    t_eval = linspace(t_span(1), t_span(2), 1000);
    f_raw = calc_field(t_eval, C_raw, w, t0, T_pulse);
    max_f = max(abs(f_raw));
    
    if max_f == 0
        cost = 1; return;
    end
    
    C_norm = C_raw * (max_f_bar / max_f);
    
    % Simulate OBE
    [~, rho] = ode45(@(t, y) obe_solver(t, y, C_norm, w, t0, T_pulse, Omega, T1, T2), t_span, [0; 0; -1]);
    
    % The cost is the negative of the final excited state fraction (to maximize it)
    cost = -rho(end, 3);
end

% Calculates the time-domain electric field
function f_t = calc_field(t, C, w, t0, T_pulse)
    s_t = zeros(size(t));
    for k = 1:length(w)
        s_t = s_t + real(C(k) .* exp(1i * w(k) * t));
    end
    f_t = s_t .* exp(-((t - t0)./T_pulse).^2);
end

% Evaluates the Optical Bloch Equations (OBE)
function dy = obe_solver(t, y, C, w, t0, T_pulse, Omega, T1, T2)
    % y(1) = rho_1, y(2) = rho_2, y(3) = rho_3
    f_t = calc_field(t, C, w, t0, T_pulse);
    Rabi = Omega * f_t;
    
    dy = zeros(3, 1);
    dy(1) = y(2) - y(1)/T2;
    dy(2) = -y(1) - y(2)/T2 + 2 * Rabi * y(3);
    dy(3) = -2 * Rabi * y(2) - (y(3) + 1)/T1;
end