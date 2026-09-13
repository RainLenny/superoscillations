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

[SO_signal_Baranov, angular_freqs_SO_Baranov] = generate_SO_from_dat_file("SO_baranov", 1, signal_scaling);

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {});

% --- Add Signal 1 ---
sig_configs(1).name = '\textbf{SO Baranov}';
sig_configs(1).data = SO_signal_Baranov;
sig_configs(1).freqs = angular_freqs_SO_Baranov;

% --- Add Cosine signals for each frequency in the SO signal ---
for k = 1:length(angular_freqs_SO_Baranov)
    freq = angular_freqs_SO_Baranov(k);
    [cos_sig, cos_freq] = generate_Cos_reference(freq);
    
    idx = length(sig_configs) + 1;
    sig_configs(idx).name = sprintf('\\boldmath$\\mathbf{%g}\\omega_0$', freq);
    sig_configs(idx).data = cos_sig;
    sig_configs(idx).freqs = cos_freq;
end

% Normalize both signals symbolically by the peak of the first signal
sig_configs = normalize_sig_configs(sig_configs, 'peak');

% Apply defaults (auto-colors)
sig_configs = prepare_signal_config(sig_configs);

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

for i = 1:length(sig_configs)
    disp(['Processing ' sig_configs(i).name '...']);

    % --- Run Solver (Pass t_common directly to avoid pchip interpolation) ---
    [~, rho_out1] = optical_bloch(t_common, rho_init, nu0_TLS1, params, sig_configs(i).data);
    [~, rho_out2] = optical_bloch(t_common, rho_init, nu0_TLS2, params, sig_configs(i).data);

    sig_configs(i).rho_tls1 = rho_out1; 
    sig_configs(i).rho_tls2 = rho_out2; 

    % --- Calculate J for each component (rho1, rho2, rho3) ---
    J_vals = zeros(length(t_axis), 3); 

    for comp = 1:3
        r1 = sig_configs(i).rho_tls1(:, comp);
        r2 = sig_configs(i).rho_tls2(:, comp);

        numerator_integrand   = (r1 - r2).^2;
        denominator_integrand = 0.5 * (r1.^2 + r2.^2);

        % cumtrapz is now highly accurate due to the dense 0.01 step size
        num_int = cumtrapz(t_common(idx_start:end), numerator_integrand(idx_start:end));
        den_int = cumtrapz(t_common(idx_start:end), denominator_integrand(idx_start:end));

        den_int(den_int == 0) = eps;

        J_vals(:, comp) = num_int ./ den_int;
    end

    sig_configs(i).J = J_vals; 
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

    for i = 1:length(sig_configs)
        plot(t_axis, sig_configs(i).J(:, comp), 'Color', sig_configs(i).color, 'LineWidth', lw, ...
             'DisplayName', sig_configs(i).name);
    end

    title(['\textbf{Distinguishability J for } \boldmath$' component_labels{comp} '$']);
    xlabel('\textbf{Integration Limit } \boldmath$T_2$');
    ylabel('\textbf{J Parameter}');
    grid on;
    legend('Location', 'best');
    xlim([T1_integ t_span(2)]);
    if exist('PlotUtils', 'class')
        PlotUtils.styleAxes(gca);
    end
end

% --- Plot rho_i in separate figures (3 figures) ---
% Shows BOTH TLS1 and TLS2 for each signal (solid = TLS1, dashed = TLS2)
for comp = 1:3
    figure('Color', 'w', 'Name', ['rho for ' component_labels{comp}]);
    hold on;

    for i = 1:length(sig_configs)
        r_tls1 = sig_configs(i).rho_tls1(:, comp);
        r_tls2 = sig_configs(i).rho_tls2(:, comp);

        % Solid line for TLS1, matching signal color
        plot(t_common, r_tls1, 'Color', sig_configs(i).color, 'LineWidth', lw, ...
            'DisplayName', sprintf('%s \\textbf{(TLS1)}', sig_configs(i).name));
        
        % Dashed line for TLS2, matching signal color
        plot(t_common, r_tls2, '--', 'Color', sig_configs(i).color, 'LineWidth', lw, ...
            'DisplayName', sprintf('%s \\textbf{(TLS2)}', sig_configs(i).name));
    end

    title(['\textbf{Time-domain trajectories for } \boldmath$' component_labels{comp} '$']);
    xlabel('\textbf{Time } \boldmath$t$');
    ylabel(['\textbf{Component } \boldmath$' component_labels{comp} '$']);
    grid on;
    legend('Location', 'best');
    xlim(t_span);
    if exist('PlotUtils', 'class')
        PlotUtils.styleAxes(gca);
    end
end

% --- Plot all input signals (Time Domain) ---
figure('Color', 'w', 'Name', 'Input Signals (Time)');
hold on;

for i = 1:length(sig_configs)
    f_vals = arrayfun(sig_configs(i).data, t_common);
    plot(t_common, f_vals, 'Color', sig_configs(i).color, 'LineWidth', lw, ...
        'DisplayName', sig_configs(i).name);
end

xlabel('\textbf{Time } \boldmath$[2\pi/\omega_0]$');
ylabel('\textbf{Amplitude [arb]}');
title('\textbf{Time-domain comparison of all input signals}');
grid on;
legend('Location', 'best');
xlim(t_span);
if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end

% --- Plot all input signals (Frequency Domain / FFT) ---
% SAMPLING Constants
f_sampling = 60/(2*pi);
% Calculate fundamental period of all aggregated frequencies
all_freqs = [];
for i = 1:length(sig_configs)
    all_freqs = [all_freqs, sig_configs(i).freqs];
end
T_period = compute_fundamental_period(all_freqs, f_sampling);
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

for i = 1:length(sig_configs)
    % Evaluate signal function over the dedicated FFT time axis
    sampled_signal = arrayfun(sig_configs(i).data, t_axis_fft);
    
    % Compute FFT
    fft_mag = fftshift(abs(fft(sampled_signal, N_fft))) / N_fft;
    
    plot(freq_axis, fft_mag, '-', 'Color', sig_configs(i).color, 'LineWidth', lw, ...
        'DisplayName', sig_configs(i).name);
end

xlabel('\textbf{Angular frequency } \boldmath$[\omega_0]$');
ylabel('\textbf{Amplitude [arb]}');
title('\textbf{Frequency-domain comparison of all input signals}');
grid on;
legend('Location', 'best');
if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(gca);
end
