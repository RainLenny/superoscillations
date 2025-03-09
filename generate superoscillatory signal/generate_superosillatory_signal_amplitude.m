% Define the sampled times t
t_span = 14;
T = linspace(-t_span, t_span, 1000); % Adjust time interval and sampling points as needed

% Define the frequencies omega_n
omega_n = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65];

% Define the target sinc function to approximate
target_function = sinc(T/pi); % sinc(x) = sin(pi*x)/(pi*x) in MATLAB

% Construct the matrix M (e^(i * omega_n * t) for all sampled times t)
M = exp(1j * omega_n' * T);

% Solve for the amplitude coefficients A using the Moore-Penrose pseudo-inverse
% Use the conjugate transpose (') for complex matrices!
input_amps = pinv(M') * target_function.';

% Generate the superoscillatory function
superoscillatory_signal = real(M' * input_amps);

% Plot the results
figure;
plot(T, target_function, 'b--', 'DisplayName', 'Target function');
hold on;
plot(T, superoscillatory_signal, 'r', 'DisplayName', 'Superoscillatory signal');
hold off;
legend;
xlabel('Time (t)');
ylabel('Amplitude');
title('Superoscillatory Signal Approximation of a Sinc Function');
grid on;