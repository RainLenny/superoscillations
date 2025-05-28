%% Plot Frequency and Impulse Response via Laplace Transform
clear; clc; close all;

%% Filter Parameters
omega0  = 1;          % Corner frequency (rad/s)
filter_order       = 7;          % Filter order
epsilon = 0.04475;    % Ripple factor

%% Symbolic Setup
syms s omega t real

% Get the symbolic transfer function H(s)
H_sym = chebychev_filter_laplace(omega0, filter_order, epsilon);

% Substitute s = j*omega
H_jw = simplify( subs(H_sym, s, 1j*omega) );
% Convert to numeric function handle
H_func = matlabFunction(H_jw, 'Vars', omega);

% Impulse Response h(t) via Inverse Laplace
h_sym_time = ilaplace(H_sym, s, t) * heaviside(t);   % h(t) in closed‐form (usually includes heaviside(t))
h_time_func = matlabFunction(h_sym_time, 'Vars', t);

%% Plot magnitude and phase
w_min = 0; w_max = 5;   % adjust as needed
figure;
subplot(2,1,1);
fplot(@(w) abs(H_func(w)), [w_min, w_max], 'LineWidth', 2);
title('Magnitude Response |H(j\omega)|');
xlabel('Frequency \omega (rad/s)');
ylabel('|H(j\omega)|');
grid on;

subplot(2,1,2);
fplot(@(w) rad2deg(angle(H_func(w))), [w_min, w_max], 'LineWidth', 2);
title('Phase Response \angle H(j\omega)');
xlabel('Frequency \omega (rad/s)');
ylabel('Phase (degrees)');
grid on;


%% Plot impulse‐response magnitude and phase
t_min = -100; t_max = 100;   % choose a window t>=0
figure;
subplot(2,1,1);
fplot(@(t) abs(h_time_func(t)), [t_min, t_max], 'LineWidth', 2);
title('Impulse Response Magnitude |h(t)|');
xlabel('Time t (s)');
ylabel('|h(t)|');
grid on;

subplot(2,1,2);
fplot(@(t) rad2deg(angle(h_time_func(t))), [t_min, t_max], 'LineWidth', 2);
title('Impulse Response Phase \angle h(t)');
xlabel('Time t (s)');
ylabel('Phase (degrees)');
grid on;
