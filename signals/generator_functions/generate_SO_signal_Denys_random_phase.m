clc;
clear;

% MATLAB Code to generate a non-superoscillating Denys signal using random phases
% Defined to have the same constituent frequencies and amplitudes (magnitudes)
% but without the phase alignment required for destructive interference.

%% 1. Define the frequencies and parameters
% Define the angular frequencies directly (rad/s)
angular_freqs_SO = [0.60, 0.70, 0.80, 0.90]; 
N = length(angular_freqs_SO);

%% 2. Generate Random Time Delays (Phases)
% Use a fixed seed for reproducible random phases
rng(1);
tau_rand = rand(1, N); 

% Define complex amplitudes
amps_SO = exp(-1i * angular_freqs_SO .* tau_rand);

% Format the variables for saving (matching other scripts)
amps_SO = amps_SO(:).';
angular_freqs_SO = angular_freqs_SO(:)';

%% 3. Save Parameters Directly
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(fullfile(output_dir, 'SO_Denys_random_phase.mat'), 'amps_SO', 'angular_freqs_SO');
disp('Parameters saved to SO_Denys_random_phase.mat');

%% 4. Logic and Plotting (Reconstruction & Visualization)
% Define a wider time range
t_continuous = linspace(-100, 100, 50001);

E_total = zeros(size(t_continuous));

% Highest frequency component for comparison (0.90 THz)
E_highest_freq = cos(max(angular_freqs_SO) * (t_continuous - tau_rand(end)));

% Reconstruct the signal: s(t) = Re( sum( c_n * exp(i * omega_n * t) ) )
for n = 1:N
    E_total = E_total + real(amps_SO(n) * exp(1i * angular_freqs_SO(n) * t_continuous));
end

% Plot 1: Combined Waveform
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1);
plot(t_continuous, E_total, 'b-', 'LineWidth', 1.5);
hold on;
title('Combined Waveform $E(t)$ (No SO)', 'Interpreter', 'latex');
xlabel('Time');
ylabel('Amplitude (arb. u.)');
grid on;
legend('Total Field', 'Location', 'best');

% Subplot 2: Zoom in on the central window
subplot(2,1,2);
plot(t_continuous, E_total, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_continuous, E_highest_freq, 'g-', 'LineWidth', 1);
title('Signal vs Highest Frequency Component ($0.9$ THz)', 'Interpreter', 'latex');
xlabel('Time ');
ylabel('Amplitude (arb. u.)');
grid on;
legend('Non-Superoscillating Field', '0.9 THz Component', 'Location', 'best');

% Display the calculated random time delays and coefficients
disp('Random Phase Assignment Complete. The time delays (tau_j) are:');
for i = 1:N
    fprintf('tau_%d (omega = %.2f rad/s) = %7.4f s\n', i, angular_freqs_SO(i), tau_rand(i));
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
title('SO Denys (No SO): Wave and Analytical Instantaneous Frequency');
xlabel('Time');
grid on;
legend('Location', 'northeast');

if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end