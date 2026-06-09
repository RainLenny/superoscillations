%% Setup Parameters
clear; clc; close all;

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

%% 3. Generate Signals & Build Struct Array
signal_scaling = 7;

% Assuming these helper functions are defined elsewhere in your path
[SO_signal_Denys, angular_freqs_SO_Denys] = generate_SO_from_dat_file('SO_Denys', 1, signal_scaling, false);
[SO_signal_Denys_no_SO, ~] = generate_SO_from_dat_file('SO_Denys_no_SO', 1, signal_scaling, false);
[SO_signal_flat, angular_freqs_SO_flat] = generate_equal_spread(angular_freqs_SO_Denys, 1, signal_scaling, false);

[Cos_signal_9, ~] = generate_Cos_reference(0.9, 1, false);
[Cos_signal_8, ~] = generate_Cos_reference(0.8, 1, false);
[Cos_signal_7, ~] = generate_Cos_reference(0.7, 1, false);
[Cos_signal_6, ~] = generate_Cos_reference(0.6, 1, false);

% Normalize both signals symbolically by the peak of the first signal
norm_sigs = normalize_signals({SO_signal_Denys, SO_signal_Denys_no_SO, SO_signal_flat, Cos_signal_9, Cos_signal_8, Cos_signal_7, Cos_signal_6}, 'peak');

% Initialize the signals struct array
signals = struct('name', {}, 'data', {}, 'y_filt1', {}, 'y_filt2', {}, 'J', {});

% --- Add Signals ---
signals(1).name = '\textbf{SO Denys}';
signals(1).data = norm_sigs{1};

signals(2).name = '\textbf{SO Denys (No SO)}';
signals(2).data = norm_sigs{2};

signals(3).name = '\textbf{flat spectrum}';
signals(3).data = norm_sigs{3};

signals(4).name = '\boldmath$\mathrm{0.9}$';
signals(4).data = norm_sigs{4};

signals(5).name = '\boldmath$\mathrm{0.8}$';
signals(5).data = norm_sigs{5};

signals(6).name = '\boldmath$\mathrm{0.7}$';
signals(6).data = norm_sigs{6};

signals(7).name = '\boldmath$\mathrm{0.6}$';
signals(7).data = norm_sigs{7};

%% 4. Main Processing Loop
dt_common = 0.01;
t_common = (t_span(1):dt_common:t_span(2))';

% Find the index corresponding to the center point
idx_center = find(t_common >= T_center, 1);
if isempty(idx_center)
    error('T_center is outside the defined time span.');
end

% --- UPDATED: Limit maximum symmetric expansion ---
% 1. Determine bounds based on array edges
max_k_bounds = min(idx_center - 1, length(t_common) - idx_center);

% 2. Determine bounds based on user-defined max window length
max_k_window = floor((max_window_length / 2) / dt_common);

% 3. Apply the strictest limit
max_k = min(max_k_bounds, max_k_window);

% Array of total window sizes (expanding evenly to left and right)
window_sizes = 2 * (0:max_k)' * dt_common;

% Indices representing the left and right boundaries of the growing window
left_indices = idx_center - (0:max_k)';
right_indices = idx_center + (0:max_k)';

for i = 1:length(signals)
    disp(['Processing ' signals(i).name '...']);

    % Evaluate input signal over common time axis
    u_in = arrayfun(signals(i).data, t_common);

    % --- Run Filter Simulation ---
    % 1. Ideal filter output
    y1 = lsim(sys, u_in, t_common);
    
    % 2. "Noisy" filter output (Add simulated noise to the output)
    y2 = y1 + time_noise_std * randn(size(y1));

    signals(i).y_filt1 = y1; 
    signals(i).y_filt2 = y2; 

    % --- Calculate J for the filter outputs ---
    numerator_integrand   = (y1 - y2).^2;
    denominator_integrand = 0.5 * (y1.^2 + y2.^2);

    % Precalculate cumulative integrals across the whole time span
    num_cum = cumtrapz(t_common, numerator_integrand);
    den_cum = cumtrapz(t_common, denominator_integrand);

    % Evaluate the integral strictly within the expanding boundaries
    num_int = num_cum(right_indices) - num_cum(left_indices);
    den_int = den_cum(right_indices) - den_cum(left_indices);

    % Prevent division by zero
    den_int(den_int == 0) = eps;

    signals(i).J = num_int ./ den_int; 
end

%% 5. Plotting
if exist('PlotUtils', 'class')
    PlotUtils.setupDefaults();
end

% Standardize line width across all plots to match Bloch script
lw = 5.0;

% --- Plot J separately ---
figure('Color', 'w', 'Name', 'Distinguishability J');
hold on;
for i = 1:length(signals)
    plot(window_sizes, signals(i).J, 'LineWidth', lw, ...
         'DisplayName', signals(i).name);
end
title('Distinguishability J (Ideal vs Noisy Filter)');
xlabel('Observation Window');
ylabel('J Parameter');
grid on;
legend('Location', 'best', 'Interpreter', 'latex');
xlim([0 max_window_length]);


% --- Plot filter outputs in a separate figure ---
figure('Color', 'w', 'Name', 'Filter Outputs');
hold on;
for i = 1:length(signals)
    % Solid line for Ideal Filter
    plot(t_common, signals(i).y_filt1, 'LineWidth', lw, ...
        'DisplayName', sprintf('%s (Ideal)', signals(i).name));
    
    % Dashed line for Noisy Filter
    plot(t_common, signals(i).y_filt2, '--', 'LineWidth', lw, ...
        'DisplayName', sprintf('%s (Noisy)', signals(i).name));
end
title('Time-domain trajectories: Ideal vs Noisy Filter');
xlabel('Time t');
ylabel('Output Amplitude');
grid on;
legend('Location', 'best', 'Interpreter', 'latex');
xlim(t_span);


% --- Plot all input signals (Time Domain) ---
figure('Color', 'w', 'Name', 'Input Signals (Time)');
hold on;
for i = 1:length(signals)
    f_vals = arrayfun(signals(i).data, t_common);
    plot(t_common, f_vals, 'LineWidth', lw, ...
        'DisplayName', signals(i).name);
end
xlabel('\boldmath$\mathrm{Time \ [2\pi/\omega_0]}$', 'Interpreter', 'latex');
ylabel('\boldmath$\mathrm{Amplitude \ [arb]}$', 'Interpreter', 'latex');
title('Time-domain comparison of all input signals');
grid on;
legend('Location', 'best', 'Interpreter', 'latex');
xlim(t_span);
if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end


% --- Plot all input signals (Frequency Domain / FFT) ---
f_sampling = 60/(2*pi);

% Calculate fundamental period of all aggregated frequencies
T_period = compute_fundamental_period([angular_freqs_SO_Denys], f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample-interval in seconds

% Setup precise FFT time axis
t_axis_fft = -signals_duration/2 : dt : signals_duration/2;
t_axis_fft = t_axis_fft(1:end-1);
t_axis_fft = t_axis_fft(:);

N_fft = length(t_axis_fft);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N_fft) * 2 * pi; % Frequency axis

figure('Color', 'w', 'Name', 'Input Signals & Filters (FFT)');
hold on;

% 1. Plot the FFT magnitudes of the input signals
for i = 1:length(signals)
    sampled_signal = arrayfun(signals(i).data, t_axis_fft);
    fft_mag = fftshift(abs(fft(sampled_signal, N_fft))) / N_fft;
    
    plot(freq_axis, fft_mag, '-', 'LineWidth', lw, ...
        'DisplayName', signals(i).name);
end

% 2. Calculate and plot the Filter Frequency Responses
H1_mag = abs(freqs(num, den, freq_axis));

% --- ADD NOISE TO THE FREQUENCY RESPONSE PLOT ---
% We add noise directly to the magnitude curve for plotting purposes
H2_mag = H1_mag + time_noise_std * randn(size(H1_mag));
H2_mag(H2_mag < 0) = 0; % Prevent the noisy magnitude from dropping below zero physically

% Plot the filter shapes using distinct dark, dashed/dotted lines 
plot(freq_axis, H1_mag, '--k', 'LineWidth', 5.0, ...
    'DisplayName', 'Ideal Filter');
plot(freq_axis, H2_mag, ':k', 'LineWidth', 2.0, ...  % Lowered line width slightly so noise is visible
    'DisplayName', 'Noisy Filter');

xlabel('\boldmath$\mathrm{Angular \ frequency \ [\omega_0]}$', 'Interpreter', 'latex');
ylabel('\boldmath$\mathrm{Amplitude \ [arb]}$', 'Interpreter', 'latex');
title('Frequency-domain: Input signals vs. Ideal and Noisy Filter Roll-off');
grid on;
legend('Location', 'best', 'Interpreter', 'latex');

% Zoom in on the x-axis to clearly see the transition band 
xlim([-2 2]);

% --- Plot Wave and Instantaneous Frequency (d_angle/dt) ---
figure('Color', 'w', 'Name', 'Instantaneous Frequency Analysis');

% Extract the time-domain data for the first three signals 
y_denys   = arrayfun(signals(1).data, t_common);
y_denys_no_SO = arrayfun(signals(2).data, t_common);
y_flat    = arrayfun(signals(3).data, t_common);

% Compute instantaneous frequency for SO Denys
z_denys = hilbert(y_denys);                       % Analytic signal
phase_denys = unwrap(angle(z_denys));             % Continuous phase
inst_freq_denys = gradient(phase_denys, dt_common); % d(angle)/dt

% Compute instantaneous frequency for SO Denys (no SO)
z_denys_no_SO = hilbert(y_denys_no_SO);
phase_denys_no_SO = unwrap(angle(z_denys_no_SO));
inst_freq_denys_no_SO = gradient(phase_denys_no_SO, dt_common);

% Compute instantaneous frequency for SO Flat
z_flat = hilbert(y_flat);                         % Analytic signal
phase_flat = unwrap(angle(z_flat));               % Continuous phase
inst_freq_flat = gradient(phase_flat, dt_common); % d(angle)/dt

% Create visual limits similar to the provided reference image
x_limits = [-5 5];
y_limits = [-4 5];

% Subplot 1: SO Denys
subplot(3,1,1);
hold on;
plot(t_common, y_denys, 'LineWidth', 2, 'DisplayName', 'wave');
plot(t_common, inst_freq_denys, 'LineWidth', 2, 'DisplayName', 'd\_angle/dt');
title('SO Denys: Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end

% Subplot 2: SO Denys (No SO)
subplot(3,1,2);
hold on;
plot(t_common, y_denys_no_SO, 'LineWidth', 2, 'DisplayName', 'wave');
plot(t_common, inst_freq_denys_no_SO, 'LineWidth', 2, 'DisplayName', 'd\_angle/dt');
title('SO Denys (No SO): Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end

% Subplot 3: Flat Spectrum
subplot(3,1,3);
hold on;
plot(t_common, y_flat, 'LineWidth', 2, 'DisplayName', 'wave');
plot(t_common, inst_freq_flat, 'LineWidth', 2, 'DisplayName', 'd\_angle/dt');
title('Flat Spectrum: Wave and Instantaneous Frequency');
xlabel('time');
xlim(x_limits);
ylim(y_limits);
grid on;
legend('Location', 'northeast');
if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end