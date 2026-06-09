%% Setup Parameters
clear; clc; close all;

rng(3);
%% 1. Global constants and Time
t_span = [-50 50];

% Integration window parameters for J calculation
T_center = 0;   % CHOSEN CENTER POINT: The point from which the window expands evenly

% For J integral calculation 
max_window_length = 7; % Total window size in time units (e.g., seconds)

%% 2. Filter Design & Noise Parameters
% We design a single analog Chebyshev Type I low-pass filter
filter_order = 4;
ripple_dB = 0.5; % Peak-to-peak passband ripple in dB

% The transition frequency is set around 0.7 to 1.0 rad/s
wc = 1.1; % Cutoff frequency for the Filter (rad/s)

% Generate continuous-time transfer functions
[num, den] = cheby1(filter_order, ripple_dB, wc, 's');
sys = tf(num, den);

% --- NOISE PARAMETERS ---
% Tune these to control how "noisy" the second filter appears
time_noise_std = 0.05;  % Standard deviation of noise added to the time-domain output

%% 3. Run Core Simulation and Plotting
params = struct();
params.time_noise_std = time_noise_std;
params.labels.filter1 = 'Ideal';
params.labels.filter2 = 'Noisy';
params.titles.distinguishability = 'Distinguishability J (Ideal vs Noisy Filter)';
params.titles.outputs = 'Time-domain trajectories: Ideal vs Noisy Filter';
params.titles.fft = 'Frequency-domain: Input signals vs. Ideal and Noisy Filter Roll-off';
params.duration_factor = 5;
params.t_span = t_span;
params.T_center = T_center;
params.max_window_length = max_window_length;
params.signal_scaling = 7;

run_filter_analysis_core(sys, [], params);