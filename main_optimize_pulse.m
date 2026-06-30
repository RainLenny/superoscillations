%% Script to optimize signal amplitudes for maximal excitation
clear; clc; close all;
addpath(genpath(pwd)); % Ensure all subfolders are on the path




NOTES:
THE NEW SIGNAL IS SUPEROSCILLATING
THE EXCITAITON DOES NOT SEEM MUCH DIFFERENT
NEED TO MAKE THE OPTIMIZED SIGNAL MUCH BETTER 
NEED TO OPTIMIZE ON VASTLY DIFFERENT DRIVE VALUES 

%% 1. Configuration & Original Signal Parameters
normalization_type = 'energy'; % Options: 'peak' or 'energy'
signal_scaling = 7;
% We generate the original SO signal to get its frequencies, amplitudes, and peak amplitude
[SO_signal_orig, angular_freqs, amps_orig] = generate_SO_from_dat_file("SO_baranov", 1, signal_scaling);

% Find the peak and energy of the original signal to use as a constraint
t_samples = linspace(0, 600, 2000);
orig_sig_vals = arrayfun(SO_signal_orig, t_samples);

max_peak_orig = max(abs(orig_sig_vals));
energy_orig = trapz(t_samples, abs(orig_sig_vals).^2);

if strcmp(normalization_type, 'peak')
    fprintf('Constraint: Peak amplitude <= %.4f\n', max_peak_orig);
elseif strcmp(normalization_type, 'energy')
    fprintf('Constraint: Energy <= %.4f\n', energy_orig);
else
    error('Unknown normalization_type: %s', normalization_type);
end

%% 2. System Parameters (from main_Bloch.m)
params.T1         = 500;
params.T2         = 300;
params.rho30      = -1.0;
params.OmegaTilde = 0.01;
nu0_target        = 1.00;

%% 3. Setup Optimization
% Optimization Variables: x = [real(amps); imag(amps)]
N_freqs = length(angular_freqs);
x0 = [real(amps_orig), imag(amps_orig)];

% Define objective function: we want to MINIMIZE -rho_3(end), which MAXIMIZES rho_3(end)
obj_fun = @(x) cost_function(x, angular_freqs, params, nu0_target);

% Define non-linear constraint based on chosen normalization
if strcmp(normalization_type, 'peak')
    nonlcon = @(x) peak_constraint(x, angular_freqs, max_peak_orig);
else
    nonlcon = @(x) energy_constraint(x, angular_freqs, energy_orig);
end

% Options for fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 300, ...
    'UseParallel', false); % Set to true if you have Parallel Computing Toolbox

%% 4. Run Optimization
fprintf('\nStarting optimization of %d complex amplitudes...\n', N_freqs);
tic;
[x_opt, fval, exitflag, output] = fmincon(obj_fun, x0, [], [], [], [], [], [], nonlcon, options);
toc;

%% 5. Results & Validation
A_opt = x_opt(1:N_freqs) + 1i * x_opt(N_freqs+1:end);

% Create optimized signal function (matching logic in generate_signal_base.m)
opt_sig_func = @(t) real(sum(conj(A_opt(:).') .* exp(1i * angular_freqs(:).' .* t(:)), 2)) .* exp(-(t(:)-250).^2./100^2);

% Simulate original response
fprintf('\nSimulating final responses...\n');
[t_orig, rho_orig] = optical_bloch([0 600], [0; 0; -1], nu0_target, params, SO_signal_orig);

% Simulate optimized response
[t_opt, rho_opt] = optical_bloch([0 600], [0; 0; -1], nu0_target, params, opt_sig_func);

fprintf('\n=== Results ===\n');
fprintf('Original Excitation (rho_3 at end): %.6f\n', rho_orig(end, 3));
fprintf('Optimized Excitation (rho_3 at end): %.6f\n', rho_opt(end, 3));

%% 6. Plotting
PlotUtils.setupDefaults();

% Time-domain Signal Comparison
figure('Color', 'w', 'Name', 'Signal Comparison');
hold on;
plot(t_samples, orig_sig_vals, 'b', 'LineWidth', 3, 'DisplayName', 'Original SO Signal');
plot(t_samples, arrayfun(opt_sig_func, t_samples), 'r', 'LineWidth', 3, 'DisplayName', 'Optimized Signal');
xlabel('Time [arb]');
ylabel('Amplitude [arb]');
title('Pulse Shape Comparison');
legend('Location', 'best');
grid on;

% Frequency-domain Amplitude Comparison
figure('Color', 'w', 'Name', 'Frequency Amplitudes');
hold on;
stem(angular_freqs, abs(amps_orig), 'b', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Original |A|');
stem(angular_freqs, abs(A_opt), 'r', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Optimized |A|');
xlabel('Angular Frequency \omega');
ylabel('Absolute Amplitude |A|');
title('Frequency Component Amplitudes');
legend('Location', 'best');
grid on;

% Excitation Trajectory
figure('Color', 'w', 'Name', 'Excitation Trajectory');
hold on;
plot(t_orig, rho_orig(:, 3), 'b', 'LineWidth', 3, 'DisplayName', 'Original \rho_3');
plot(t_opt, rho_opt(:, 3), 'r', 'LineWidth', 3, 'DisplayName', 'Optimized \rho_3');
xlabel('Time [arb]');
ylabel('\rho_3 (Excitation)');
title('Target TLS Excitation Trajectory');
legend('Location', 'best');
grid on;

% --- Instantaneous Frequency Analysis ---
sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'tgrid', {}, 'Pe', {});

sig_configs(1).name = 'Original SO Signal';
sig_configs(1).data = SO_signal_orig;
sig_configs(1).color = 'b';
sig_configs(1).tgrid = t_orig;
sig_configs(1).Pe = (1 + rho_orig(:, 3)) / 2;

sig_configs(2).name = 'Optimized Signal';
sig_configs(2).data = opt_sig_func;
sig_configs(2).color = 'r';
sig_configs(2).tgrid = t_opt;
sig_configs(2).Pe = (1 + rho_opt(:, 3)) / 2;

T_final = 600;
t_inst = linspace(0, T_final, 30000)'; 
dt_inst = t_inst(2) - t_inst(1);

for i = 1:length(sig_configs)
    y_sig = arrayfun(sig_configs(i).data, t_inst);
    inst_freq = compute_instantaneous_frequency(y_sig, dt_inst);
    Pe_inst = interp1(sig_configs(i).tgrid, sig_configs(i).Pe, t_inst);
    
    figure('Color', 'w', 'Name', sprintf('%s: Full Analysis', sig_configs(i).name));
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    
    ax(1) = nexttile; 
    plot(t_inst, Pe_inst, 'LineWidth', 2.5, 'Color', sig_configs(i).color);
    ylabel('P_e'); title(sprintf('%s: Excitation & Instantaneous Frequency', sig_configs(i).name));
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


%% Helper Functions

function f = cost_function(x, freqs, params, nu0)
    N = length(freqs);
    A = x(1:N) + 1i * x(N+1:end);
    
    % Fast inline signal generation for ODE scalar 't'
    sig_func = @(t) real(sum(conj(A(:).') .* exp(1i * freqs(:).' .* t), 2)) .* exp(-(t-250).^2./100^2);
    
    rho_init = [0; 0; -1];
    t_span = [0 600];
    
    % Solve Bloch equations
    [~, rho] = optical_bloch(t_span, rho_init, nu0, params, sig_func);
    
    % Maximize rho_3 at the end of the simulation
    f = -rho(end, 3);
end

function [c, ceq] = peak_constraint(x, freqs, max_peak)
    N = length(freqs);
    A = x(1:N) + 1i * x(N+1:end);
    
    % Sample time window around the Gaussian peak (250 +/- 150 is where the peak should be)
    t = linspace(100, 400, 400)'; 
    sig = real(sum(conj(A(:).') .* exp(1i * freqs(:).' .* t), 2)) .* exp(-(t-250).^2./100^2);
    
    % Constraint: max(abs(sig)) - max_peak <= 0
    c = max(abs(sig)) - max_peak;
    ceq = [];
end

function [c, ceq] = energy_constraint(x, freqs, max_energy)
    N = length(freqs);
    A = x(1:N) + 1i * x(N+1:end);
    
    % Sample time window around the Gaussian peak (250 +/- 150)
    t = linspace(100, 400, 400)'; 
    sig = real(sum(conj(A(:).') .* exp(1i * freqs(:).' .* t), 2)) .* exp(-(t-250).^2./100^2);
    
    % Constraint: energy - max_energy <= 0
    energy = trapz(t, abs(sig).^2);
    c = energy - max_energy;
    ceq = [];
end
