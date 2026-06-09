%% Setup Parameters
clear; clc; close all;

%% 1. Global constants and Time
t_span = [-50 50];

% Integration window parameters for J calculation
T_center = 0;   % CHOSEN CENTER POINT: The point from which the window expands evenly

%For J integral calculation 
max_window_length =7; % Total window size in time units (e.g., seconds)

%% 2. Filter Design (Replaces TLS Parameters)
% We design two very similar analog Chebyshev Type I low-pass filters
filter_order = 4;
ripple_dB = 0.5; % Peak-to-peak passband ripple in dB

% % The transition frequency is set around 0.7 to 1.0 rad/s
% wc1 = 2; % Cutoff frequency for Filter 1 (rad/s)
% wc2 = 2.01; % Cutoff frequency for Filter 2 (rad/s)

% The transition frequency is set around 0.7 to 1.0 rad/s
wc1 = 3; % Cutoff frequency for Filter 1 (rad/s)
wc2 = 3.05; % Cutoff frequency for Filter 2 (rad/s)

% Generate continuous-time transfer functions
[num1, den1] = cheby1(filter_order, ripple_dB, wc1, 's');
[num2, den2] = cheby1(filter_order, ripple_dB, wc2, 's');

sys1 = tf(num1, den1);
sys2 = tf(num2, den2);

%% 3. Run Core Simulation and Plotting
params = struct();
params.time_noise_std = 0;
params.labels.filter1 = 'Filter 1';
params.labels.filter2 = 'Filter 2';
params.titles.distinguishability = 'Distinguishability J for Filter Outputs';
params.titles.outputs = 'Time-domain trajectories of Filter Outputs';
params.titles.fft = 'Frequency-domain comparison of all input signals vs. Filter Roll-off';
params.duration_factor = 100;
params.t_span = t_span;
params.T_center = T_center;
params.max_window_length = max_window_length;
params.signal_scaling = 7;

run_filter_analysis_core(sys1, sys2, params);