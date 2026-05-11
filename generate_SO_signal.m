clc;
clear;

% MATLAB Code to generate a superoscillating signal
% Target: 2 periods of a cosine with frequency omega_0 = 1
% Using 5 frequencies equally spaced up to 0.7

% 1. Define the frequencies
N = 5;
omega_max = 0.7;
% Create 5 equally spaced frequencies from 0.14 to 0.70
angular_freqs_SO = linspace(omega_max/N, omega_max, N); 

% 2. Define the target constraint points
% Two periods of a cosine with omega_0 = 1 span from t = 0 to 4*pi.
% We select 5 points corresponding to the peaks and troughs:
omega_0 = 1;
t_j = (0:N-1) * (pi / omega_0); % Results in [0, pi, 2*pi, 3*pi, 4*pi]

% Target values at these points to mimic cos(omega_0 * t)
s_j = cos(omega_0 * t_j); % Results in [1, -1, 1, -1, 1]

% 3. Set up the linear system of equations: A * c = s_j
% The matrix A is constructed using A_jn = exp(i * omega_n * t_j)
A = exp(1i * t_j' * angular_freqs_SO);

% Solve the linear system for the complex coefficients c_n
amps_SO = A \ s_j';

% 4. Reconstruct the superoscillating signal over a continuous time range
% We calculate it over a wider window to see the behavior outside the SO region
t_continuous = linspace(-5, 20, 1000);
s_t = zeros(size(t_continuous));

% Reconstruct the signal: s(t) = Re( sum( c_n * exp(i * omega_n * t) ) )
for n = 1:N
    s_t = s_t + real(amps_SO(n) * exp(1i * angular_freqs_SO(n) * t_continuous));
end

% 5. Plot the results
figure;
plot(t_continuous, s_t, 'b-', 'LineWidth', 1.5);
hold on;

% Plot the target constraint points
plot(t_j, s_j, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Highlight the boundaries of the superoscillating region [0, 4*pi]
xline(0, 'k--', 'LineWidth', 1.5);
xline(4*pi, 'k--', 'LineWidth', 1.5);

title('Superoscillating Signal');
xlabel('Time');
ylabel('Amplitude');
legend('Superoscillating Signal s(t)', 'Constraint Points', 'Super Region Boundaries', 'Location', 'best');
grid on;

% Display the calculated coefficients in the command window
disp('The calculated complex coefficients c_n are:');
disp(amps_SO);


amps_SO = amps_SO(:).';
angular_freqs_SO = angular_freqs_SO(:)';
%save('SO_signal_data.mat', 'amps_SO', 'angular_freqs_SO');
%load('SO_signal_data.mat', 'amps_SO', 'angular_freqs_SO');