clc;
clear;

% MATLAB Code to generate a superoscillating signal via destructive interference
% Based on the heuristic time-delay optimization method from:
% "Superoscillations Deliver Superspectroscopy" (McCaul et al., 2023)

%% 1. Define the frequencies and parameters
% The user requested frequencies between 0.6 and 0.9. 
% Define the angular frequencies directly (rad/s)
angular_freqs_SO = [0.60, 0.70, 0.80, 0.90]; 
N = length(angular_freqs_SO);

% Define the optimization window [-T_SO, T_SO] 
% This is the region where we want destructive interference to occur.
T_SO = pi; 

% Time vector specifically for the integration/minimization window
t_opt = linspace(-T_SO, T_SO, 500000);
dt = t_opt(2) - t_opt(1);

%% 2. Define the Objective Functioz
% We want to find time delays (tau) that minimize the integrated intensity.
% I(tau) = Integral of ( sum( cos(omega_j * (t - tau_j)) ) )^2 dt
% Note: The paper sets amplitudes a_j = 1 to avoid trivial zero solutions.

% Anonymous function for the cost/objective calculation
objective = @(tau) sum( sum( cos(angular_freqs_SO' .* (t_opt - tau')), 1 ).^2 ) * dt;

%% 3. Optimize Time Delays
% The objective function is non-convex and has many local minima, meaning 
% different initial guesses will yield different superoscillatory shapes.
rng(42); % Seed for reproducibility
tau0 = rand(1, N); % Initial guess for the time delays

% Use fminsearch (Nelder-Mead simplex) for unconstrained optimization
options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000, 'Display', 'iter');
disp('Optimizing time delays to minimize intensity...');
tau_opt = fminsearch(objective, tau0, options);

% Define complex amplitudes
amps_SO = exp(-1i * angular_freqs_SO .* tau_opt);

% Format the variables for saving (matching other scripts)
amps_SO = amps_SO(:).';
angular_freqs_SO = angular_freqs_SO(:)';

%% 4. Save Parameters Directly
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'signals', 'data');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(fullfile(output_dir, 'SO_Denys.mat'), 'amps_SO', 'angular_freqs_SO');
disp('Parameters saved to SO_Denys.mat');

%% 5. Logic and Plotting (Reconstruction & Visualization)
% Define a wider time range to observe the signal inside and outside the window
t_continuous = linspace(-100, 100, 50001);

E_total = zeros(size(t_continuous));

% Highest frequency component for comparison (0.90 THz)
E_highest_freq = cos(max(angular_freqs_SO) * (t_continuous - tau_opt(end)));

for n = 1:N
    E_total = E_total + real(amps_SO(n) * exp(1i * angular_freqs_SO(n) * t_continuous));
end

% Plot 1: Combined Waveform
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1);
plot(t_continuous, E_total, 'b-', 'LineWidth', 1.5);
hold on;
% Highlight the minimization window
xline(-T_SO, 'k--', 'LineWidth', 1);
xline(T_SO, 'k--', 'LineWidth', 1);
patch([-T_SO, T_SO, T_SO, -T_SO], [-max(abs(E_total)), -max(abs(E_total)), max(abs(E_total)), max(abs(E_total))], ...
    'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
title('Combined Waveform $E(t)$', 'Interpreter', 'latex');
xlabel('Time');
ylabel('Amplitude (arb. u.)');
grid on;
legend('Total Field', 'Minimization Boundary', 'Location', 'best');

% Subplot 2: Zoom in on the superoscillatory window
subplot(2,1,2);
plot(t_continuous, E_total, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_continuous, E_highest_freq, 'g-', 'LineWidth', 1);
xline(-T_SO, 'k--', 'LineWidth', 1);
xline(T_SO, 'k--', 'LineWidth', 1);
title('Superoscillation vs Highest Frequency Component ($0.9$ THz)', 'Interpreter', 'latex');
xlabel('Time ');
ylabel('Amplitude (arb. u.)');
grid on;
legend('Superoscillating Field', '0.9 THz Component', 'Location', 'best');

% Display the calculated time delays and coefficients
disp('Optimization Complete. The calculated time delays (tau_j) are:');
for i = 1:N
    fprintf('tau_%d (omega = %.2f rad/s) = %7.4f s\n', i, angular_freqs_SO(i), tau_opt(i));
end

disp('The calculated complex coefficients c_n are:');
for i=1:N
    fprintf('c_%d = %7.3f %+.3fi\n', i, real(amps_SO(i)), imag(amps_SO(i)));
end

% Analytical Instantaneous Frequency Analysis
% Preallocate arrays for the analytic signal and its derivative
z_analytical = zeros(size(t_continuous));
dz_analytical = zeros(size(t_continuous));

% Construct the exact analytic signal z(t) and its exact derivative z'(t) using stored negated amps_SO
for n = 1:N
    term = amps_SO(n) * exp(1i * angular_freqs_SO(n) * t_continuous);
    z_analytical = z_analytical + term;
    dz_analytical = dz_analytical + 1i * angular_freqs_SO(n) * term;
end

inst_freq_analytical = imag(dz_analytical ./ z_analytical); 

figure('Color', 'w', 'Name', 'Analytical Instantaneous Frequency Analysis');
plot(t_continuous, real(z_analytical), 'LineWidth', 2, 'DisplayName', 'wave s(t)'); % real(z_analytical) is identical to s_t
hold on;
plot(t_continuous, inst_freq_analytical, 'LineWidth', 2, 'DisplayName', 'Analytical \omega_{inst}(t)');
title('SO Denys: Wave and Analytical Instantaneous Frequency');
xlabel('Time');
grid on;
legend('Location', 'northeast');

if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end