% =========================================================================
% Optimal Control of a Two-Level System in a Superoscillating Field
% Based on Baranov et al. (2014) and Morzhin & Pechen (2019) https://iopscience.iop.org/article/10.1070/RM9835/pdf?casa_token=-3XATwVwsE0AAAAA:x_KLLmNGWrJyrfOw6OlXm5bZVhG-JhTxa5_QYBkxcmMj0Iv6-gfpkW3XEhgejymP8Yl3M79JMZxS3_iTlT5VkSLItw
% =========================================================================

clear; clc; close all;

%% 1. Define Parameters
T1 = 500;           % Relaxation time
T2 = 300;           % Dephasing time
Omega = 0.01;       % Rabi constant
w0 = 1;             % TLS transition frequency (normalized)
T_pulse = 100;      % Pulse duration parameter
t0_pulse = 200;     % Pulse center
tspan = [0 600];    % Simulation time
y0 = [0; 0; -1];    % Initial state: [rho1, rho2, rho3] -> ground state (rho3 = -1)

% Frequencies used in the superoscillating signal
n_freq = 1:5;
w = 0.18 * n_freq;  % w_n = 0.18, 0.36, 0.54, 0.72, 0.90

% Original complex coefficients from Baranov et al. Eq. (6)
c_orig = [-0.156 + 0.331i, ...
          -0.862 - 1.042i, ...
           2.341 - 0.601i, ...
          -0.502 + 2.634i, ...
          -1.820 - 1.322i];

%% 2. Original Signal Setup & Normalization
A_orig = 7; % Amplitude scalar used in the paper
t_plot = linspace(tspan(1), tspan(2), 2000);

% Function to evaluate the field
eval_field = @(t, c, scale) scale * exp(-(t - t0_pulse).^2 / T_pulse^2) .* ...
    sum(real(reshape(c,[],1) .* exp(1i * reshape(w,[],1) .* t)), 1);

% Compute original field to find its peak amplitude
f_orig_unscaled = eval_field(t_plot, c_orig, A_orig);
peak_amp = max(abs(f_orig_unscaled)); 
fprintf('Original signal peak amplitude: %.4f\n', peak_amp);

%% 3. Optimal Control using Nelder-Mead (CRAB-style parameterization)
% The control variables are the real and imaginary parts of the 5 coefficients
disp('Optimizing signal coefficients to maximize excitation...');
x0 = [real(c_orig), imag(c_orig)]; 

% Objective function: maximizes max(rho_3) by minimizing -max(rho_3)
obj_fun = @(x) cost_function(x, w, peak_amp, tspan, y0, T1, T2, Omega, t0_pulse, T_pulse);

% Set optimization options (Nelder-Mead simplex method)
options = optimset('Display', 'iter', 'MaxIter', 400, 'MaxFunEvals', 1000);
x_opt = fminsearch(obj_fun, x0, options);

% Reconstruct optimized complex coefficients
c_opt = x_opt(1:5) + 1i * x_opt(6:10);

%% 4. Simulate OBE for Both Signals
disp('Simulating final dynamics for both signals...');
% Find proper scaling for the optimal signal to enforce the exact same peak amplitude
f_opt_test = eval_field(t_plot, c_opt, 1);
scale_opt = peak_amp / max(abs(f_opt_test));

% Solve OBE for original
ode_orig = @(t,y) bloch_eq(t, y, c_orig, w, A_orig, T1, T2, Omega, t0_pulse, T_pulse);
[t_orig, y_orig] = ode45(ode_orig, t_plot, y0);

% Solve OBE for optimal
ode_opt = @(t,y) bloch_eq(t, y, c_opt, w, scale_opt, T1, T2, Omega, t0_pulse, T_pulse);
[t_opt, y_opt] = ode45(ode_opt, t_plot, y0);

% Calculate fields over time for plotting
f_orig_t = eval_field(t_plot, c_orig, A_orig);
f_opt_t  = eval_field(t_plot, c_opt, scale_opt);

% Excitation probability: rho_22 = (rho_3 + 1)/2
exc_orig = (y_orig(:,3) + 1) / 2;
exc_opt  = (y_opt(:,3) + 1) / 2;

fprintf('\nMax Excitation (Original): %.2f%%\n', max(exc_orig)*100);
fprintf('Max Excitation (Optimized): %.2f%%\n', max(exc_opt)*100);

%% 5. Plotting
figure('Position', [100, 100, 1200, 400], 'Name', 'Optimal Control of Superoscillating Field');

% --- A. Time-domain Signals ---
subplot(1,3,1);
plot(t_plot, f_orig_t, 'b', 'LineWidth', 1.2); hold on;
plot(t_plot, f_opt_t, 'r--', 'LineWidth', 1.2);
title('Electric Field (Time Domain)');
xlabel('Time (\omega_0^{-1})');
ylabel('Electric field f(t)');
legend('Original (Baranov)', 'Optimized', 'Location', 'best');
grid on; axis tight;

% --- B. Frequency-domain Spectra ---
subplot(1,3,2);
Fs = 1 / (t_plot(2) - t_plot(1));
L = length(t_plot);
f_omega = (Fs * (0:(L/2)) / L) * 2 * pi; % Angular frequency axis

% FFT original
Y_orig = fft(f_orig_t);
P1_orig = abs(Y_orig/L);
P1_orig = P1_orig(1:floor(L/2)+1); P1_orig(2:end-1) = 2*P1_orig(2:end-1);

% FFT optimized
Y_opt = fft(f_opt_t);
P1_opt = abs(Y_opt/L);
P1_opt = P1_opt(1:floor(L/2)+1); P1_opt(2:end-1) = 2*P1_opt(2:end-1);

plot(f_omega, P1_orig, 'b', 'LineWidth', 1.5); hold on;
plot(f_omega, P1_opt, 'r--', 'LineWidth', 1.5);
xline(w0, 'k:', 'LineWidth', 1.5, 'Label', '\omega_0 (Resonance)');
xlim([0, 1.2]);
title('Spectral Density');
xlabel('Frequency \omega');
ylabel('Amplitude');
legend('Original', 'Optimized', 'Location', 'best');
grid on;

% --- C. Excitation Probability ---
subplot(1,3,3);
plot(t_orig, exc_orig, 'b', 'LineWidth', 1.5); hold on;
plot(t_opt, exc_opt, 'r', 'LineWidth', 1.5);
title('Excitation Probability (\rho_{22})');
xlabel('Time (\omega_0^{-1})');
ylabel('Probability');
legend('Original', 'Optimized', 'Location', 'best');
grid on; axis tight;
ylim([0 1.05]);

%% Helper Functions

function neg_exc = cost_function(x, w, target_peak, tspan, y0, T1, T2, Omega, t0_pulse, T_pulse)
    % Reconstruct complex coefficients
    c = x(1:5) + 1i * x(6:10);
    
    % Ensure peak constraint
    t_test = linspace(tspan(1), tspan(2), 500);
    f_test = exp(-(t_test - t0_pulse).^2 / T_pulse^2) .* ...
             sum(real(reshape(c,[],1) .* exp(1i * reshape(w,[],1) .* t_test)), 1);
    
    peak_current = max(abs(f_test));
    if peak_current < 1e-6
        neg_exc = 0; % Penalize dead signals
        return;
    end
    scale = target_peak / peak_current;
    
    % Simulate OBE
    odefun = @(t,y) bloch_eq(t, y, c, w, scale, T1, T2, Omega, t0_pulse, T_pulse);
    [~, Y] = ode45(odefun, tspan, y0);
    
    % We want to maximize the maximum excitation probability achieved during the pulse
    max_rho3 = max(Y(:,3));
    neg_exc = -max_rho3; 
end

function dydt = bloch_eq(t, y, c, w, scale, T1, T2, Omega, t0_pulse, T_pulse)
    % Calculate the field f(t) at time t
    s_t = sum(real(c(:) .* exp(1i * w(:) * t)));
    f_val = scale * s_t * exp(-(t - t0_pulse)^2 / T_pulse^2);
    
    rho1 = y(1);
    rho2 = y(2);
    rho3 = y(3);
    
    % Optical Bloch Equations (Assuming w0 = 1, as per Baranov eq. 5)
    dydt = zeros(3,1);
    dydt(1) = rho2 - rho1 / T2;
    dydt(2) = -rho1 - rho2 / T2 + 2 * Omega * f_val * rho3;
    dydt(3) = -2 * Omega * f_val * rho2 - (rho3 - (-1)) / T1;
end