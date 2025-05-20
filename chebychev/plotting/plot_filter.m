%% Plot Magnitude and Phase Response of a Chebyshev Filter

clear; clc; close all;

%% Filter Parameters
omega0 = 1;  % Corner frequency (rad/s)
N = 9;         % Filter order
epsilon = 0.04475;  % Ripple factor 

% Generate the Chebyshev Filter symbolic transfer function (with s = j*omega)
H_sym = chebychev_filter(omega0, N, epsilon);

% Convert the symbolic function to a numeric function handle.
% This makes it vectorized and avoids conflicts with MATLAB’s built-in 'filter'
H_func = matlabFunction(H_sym, 'Vars', sym('omega'));

%% Inverse Fourier Transform
% Here we compute the impulse response h(t) as the inverse Fourier transform of H(jω)
syms t real
syms omega
h_sym_time = ifourier(H_sym, omega, t);
DO THIS WITH LAPLASE in a different file

% Convert the symbolic impulse response to a numeric function handle for plotting.
h_time_func = matlabFunction(h_sym_time, 'Vars', t);

% Define time range for plotting the impulse response.
t_min = -100;
t_max = 100;

%% Plot the Frequency Response
% Define frequency range for plotting
w_min = 0;
w_max = 5;

figure;

% Magnitude response plot
subplot(2,1,1);
fplot(@(w) abs(H_func(w)), [w_min, w_max], 'LineWidth', 2);
% fplot(@(w) mag2db(abs(H_func(w))), [w_min, w_max], 'LineWidth', 2);
title('Magnitude Response of the Chebyshev Filter');
xlabel('Frequency (rad/s)');
ylabel('|H(j\omega)|');
grid on;

% Phase response plot (converted from radians to degrees)
subplot(2,1,2);
fplot(@(w) rad2deg(angle(H_func(w))), [w_min, w_max], 'LineWidth', 2);
title('Phase Response of the Chebyshev Filter');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');
grid on;

%% Plot Fourier
figure;
% Impulse Response Magnitude Plot
subplot(2,1,1);
fplot(@(t) abs(h_time_func(t)), [t_min, t_max], 'LineWidth', 2);
title('Impulse Response Magnitude of the Chebyshev Filter');
xlabel('Time (s)');
ylabel('|h(t)|');
grid on;

% Impulse Response Phase Plot
subplot(2,1,2);
fplot(@(t) rad2deg(angle(h_time_func(t))), [t_min, t_max], 'LineWidth', 2);
title('Impulse Response Phase of the Chebyshev Filter');
xlabel('Time (s)');
ylabel('Phase (degrees)');
grid on;
