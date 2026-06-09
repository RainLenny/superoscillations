clc;
clear;

% MATLAB Code to generate the specific superoscillating signal from the paper:
% "Abrupt Rabi oscillations in a superoscillating electric field"

% 1. Define the frequencies
N = 5;
% The paper uses five harmonic oscillations with frequencies omega_n = 0.18n
n_idx = 1:N;
angular_freqs_SO = 0.18 * n_idx;

% 2. Define the target constraint points
% Five points t_n = pi * n / omega_0, for n = 0,...,4, with omega_0 = 1
omega_0 = 1;
t_j = (0:N-1) * (pi / omega_0); % Results in [0, pi, 2*pi, 3*pi, 4*pi]

% Target values assigned to these points: s1=s3=s5=-1, s2=s4=1
s_j = [-1, 1, -1, 1, -1];

% 3. Set up the linear system of equations: A * c = s_j
% The matrix A is constructed using A_jn = exp(i * omega_n * t_j)
A = exp(1i * t_j' * angular_freqs_SO);

% Solve the linear system for the complex coefficients c_n
amps_SO = A \ s_j';

%% Save

amps_SO = amps_SO(:).';
angular_freqs_SO = angular_freqs_SO(:)';

% Save directly to the signals/data directory
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'signals', 'data');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(fullfile(output_dir, 'SO_Baranov.mat'), 'amps_SO', 'angular_freqs_SO');

%% PLOT

% 4. Reconstruct the superoscillating signal over a continuous time range
% The paper plots this from t=0 to roughly t=55 in Figure 1(b)
t_continuous = linspace(-55, 55, 2000);
s_t = zeros(size(t_continuous));

% Reconstruct the signal: s(t) = Re( sum( c_n * exp(i * omega_n * t) ) )
for n = 1:N
    s_t = s_t + real(amps_SO(n) * exp(1i * angular_freqs_SO(n) * t_continuous));
end

% 5. Plot the results to mimic Figure 1(a) and 1(b)
figure('Position', [100, 100, 800, 600]);

% Subplot for the constraint points (Figure 1a match)
subplot(2,1,1);
plot(t_j, s_j, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
yline(0, 'k-', 'Color', [0.8 0.8 0.8]);
title('Predefined pattern of points (Fig. 1a)');
xlabel('Time');
xlim([-1, 13]);
ylim([-2, 2]);

% Subplot for the SO signal (Figure 1b match)
subplot(2,1,2);
plot(t_continuous, s_t, 'b-', 'LineWidth', 1.5);
hold on;
% Plot the constraint points on the signal to verify alignment
plot(t_j, s_j, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
title('Superoscillating Function s(t) (Fig. 1b)');
xlabel('Time');
xlim([-2, 55]);
grid on;

% Display the calculated coefficients in the command window
disp('The calculated complex coefficients c_n are:');
for i=1:N
    fprintf('c_%d = %7.3f %+.3fi\n', i, real(amps_SO(i)), imag(amps_SO(i)));
end



%% 6. Analytical Instantaneous Frequency Analysis
% Preallocate arrays for the analytic signal and its derivative
z_analytical = zeros(size(t_continuous));
dz_analytical = zeros(size(t_continuous));

% Construct the exact analytic signal z(t) and its exact derivative z'(t)
for n = 1:N
    term = amps_SO(n) * exp(1i * angular_freqs_SO(n) * t_continuous);
    z_analytical = z_analytical + term;
    dz_analytical = dz_analytical + 1i * angular_freqs_SO(n) * term;
end

% REMOVED THE NEGATIVE SIGN HERE:
% Standard math identity: omega = Im(z' / z)
inst_freq_analytical = imag(dz_analytical ./ z_analytical); 

figure('Color', 'w', 'Name', 'Analytical Instantaneous Frequency Analysis');
plot(t_continuous, real(z_analytical), 'LineWidth', 2, 'DisplayName', 'wave s(t)'); 
hold on;
plot(t_continuous, inst_freq_analytical, 'LineWidth', 2, 'DisplayName', 'Analytical \omega_{inst}(t)');
title('SO Baranov: Wave and Analytical Instantaneous Frequency');
xlabel('Time');
ylim([-6 6]); 
grid on;
legend('Location', 'northeast');

if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end