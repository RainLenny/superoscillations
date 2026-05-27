%% Setup Parameters
clear; clc; close all;

%% 1. Global constants and Time
t_span = [0 600];

% Integration window parameters for J calculation
T1_integ = 235;   % fixed start time (T_start)

%% 2. Initial Conditions and System Constants
rho_init = [0; 0; -1];

params.T1         = 500;
params.T2         = 300;
params.rho30      = -1.0;
params.OmegaTilde = 0.01;

% Define the two TLS systems to compare
nu0_TLS1 = 1.00;
nu0_TLS2 = 1.01;

%% 3. Generate Signals & Build Struct Array
signal_scaling = 7;

[SO_signal_Baranov, angular_freqs_SO_Baranov] = generate_SO_Baranov(1, signal_scaling);
[SO_signal_flat, angular_freqs_SO_flat] = generate_SO_equal_spread(1, signal_scaling);
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9);

% Normalize both signals symbolically by the peak of the first signal
norm_sigs = normalize_signals({SO_signal_Baranov, SO_signal_flat, Cos_signal}, 'peak');

% Generate distinct colors for plotting

% Initialize the signals struct array
signals = struct('name', {}, 'data', {}, 'rho_tls1', {}, 'rho_tls2', {}, 'J', {});

% --- Add Signal 1 ---
signals(1).name = '\textbf{SO Baranov}';
signals(1).data    = norm_sigs{1};

% --- Add Signal 2 ---
signals(2).name = '\textbf{SO flat}';
signals(2).data    = norm_sigs{2};

% --- Add Signal 3 ---
signals(3).name = '\boldmath$\mathrm{0.9\omega_0}$';
signals(3).data    = norm_sigs{3};


%% 4. Main Processing Loop
% FIX: Ensure dt is small enough to capture the 30 rad/s high-frequency components
% avoiding both visual aliasing and integration (cumtrapz) errors.
dt_common = 0.01;
t_common = (t_span(1):dt_common:t_span(2))';

% Index corresponding to T1_integ
idx_start = find(t_common >= T1_integ, 1);
if isempty(idx_start), idx_start = 1; end

% Time axis for J plots (T2)
t_axis = t_common(idx_start:end);

for i = 1:length(signals)
    disp(['Processing ' signals(i).name '...']);

    % --- Run Solver (Pass t_common directly to avoid pchip interpolation) ---
    [~, rho_out1] = optical_bloch(t_common, rho_init, nu0_TLS1, params, signals(i).data);
    [~, rho_out2] = optical_bloch(t_common, rho_init, nu0_TLS2, params, signals(i).data);

    signals(i).rho_tls1 = rho_out1; 
    signals(i).rho_tls2 = rho_out2; 

    % --- Calculate J for each component (rho1, rho2, rho3) ---
    J_vals = zeros(length(t_axis), 3); 

    for comp = 1:3
        r1 = signals(i).rho_tls1(:, comp);
        r2 = signals(i).rho_tls2(:, comp);

        numerator_integrand   = (r1 - r2).^2;
        denominator_integrand = 0.5 * (r1.^2 + r2.^2);

        % cumtrapz is now highly accurate due to the dense 0.01 step size
        num_int = cumtrapz(t_common(idx_start:end), numerator_integrand(idx_start:end));
        den_int = cumtrapz(t_common(idx_start:end), denominator_integrand(idx_start:end));

        den_int(den_int == 0) = eps;

        J_vals(:, comp) = num_int ./ den_int;
    end

    signals(i).J = J_vals; 
end
%% 5. Plotting
component_labels = {'\rho_1', '\rho_2', '\rho_3'};
PlotUtils.setupDefaults();

% Standardize line width across all plots
lw = 5.0;

% --- Plot J_i separately (3 figures) ---
for comp = 1:3
    figure('Color', 'w', 'Name', ['J for ' component_labels{comp}]);
    hold on;

    for i = 1:length(signals)
        plot(t_axis, signals(i).J(:, comp), 'LineWidth', lw, ...
             'DisplayName', signals(i).name);
    end

    title(['Distinguishability J for ' component_labels{comp}]);
    xlabel('Integration Limit T_2');
    ylabel('J Parameter');
    grid on;
    legend('Location', 'best', 'Interpreter', 'latex');
    xlim([T1_integ t_span(2)]);
end

% --- Plot rho_i in separate figures (3 figures) ---
% Shows BOTH TLS1 and TLS2 for each signal (solid = TLS1, dashed = TLS2)
for comp = 1:3
    figure('Color', 'w', 'Name', ['rho for ' component_labels{comp}]);
    hold on;

    for i = 1:length(signals)
        r_tls1 = signals(i).rho_tls1(:, comp);
        r_tls2 = signals(i).rho_tls2(:, comp);

        % Solid line for TLS1, matching signal color
        plot(t_common, r_tls1, 'LineWidth', lw, ...
            'DisplayName', sprintf('%s (TLS1)', signals(i).name));
        
        % Dashed line for TLS2, matching signal color
        plot(t_common, r_tls2, '--', 'LineWidth', lw, ...
            'DisplayName', sprintf('%s (TLS2)', signals(i).name));
    end

    title(['Time-domain trajectories for ' component_labels{comp}]);
    xlabel('Time t');
    ylabel(['Component ' component_labels{comp}]);
    grid on;
    legend('Location', 'best', 'Interpreter', 'latex');
    xlim(t_span);
end

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
% SAMPLING Constants
f_sampling = 60/(2*pi);
% Calculate fundamental period of all aggregated frequencies
T_period = compute_fundamental_period([angular_freqs_SO_flat, angular_freqs_COS], f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample-interval in seconds

% Setup precise FFT time axis
t_axis_fft = -signals_duration/2 : dt : signals_duration/2;
t_axis_fft = t_axis_fft(1:end-1);
t_axis_fft = t_axis_fft(:);

N_fft = length(t_axis_fft);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N_fft) * 2 * pi; % Frequency axis

figure('Color', 'w', 'Name', 'Input Signals (FFT)');
hold on;

for i = 1:length(signals)
    % Evaluate signal function over the dedicated FFT time axis
    sampled_signal = arrayfun(signals(i).data, t_axis_fft);
    
    % Compute FFT
    fft_mag = fftshift(abs(fft(sampled_signal, N_fft))) / N_fft;
    
    plot(freq_axis, fft_mag, '-', 'LineWidth', lw, ...
        'DisplayName', signals(i).name);
end

xlabel('\boldmath$\mathrm{Angular \ frequency \ [\omega_0]}$', 'Interpreter', 'latex');
ylabel('\boldmath$\mathrm{Amplitude \ [arb]}$', 'Interpreter', 'latex');
title('Frequency-domain comparison of all input signals');
grid on;
legend('Location', 'best', 'Interpreter', 'latex');
