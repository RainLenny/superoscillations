%% PATH SETUP & CONSTANTS
clear; clc;
% Add all project subfolders to search path
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(project_root));

% Try to use PlotUtils if it exists, otherwise ignore
try
    PlotUtils.setupDefaults();
catch
end

%% 1. Define Signals
signal_scaling = 0.2;

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {});

% --- Signal 1: SO ---
[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, false);
sig_configs(1).name = '\textbf{SO }';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;

% --- Signal 2: COS ---
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(max(angular_freqs_SO), signal_scaling, false);
sig_configs(2).name = '\boldmath$\mathbf{0.9\omega_0}$';
sig_configs(2).data = Cos_signal;
sig_configs(2).color = 'b';
sig_configs(2).freqs = angular_freqs_COS;

% Normalize all signals symbolically by the peak of the first signal
try
    sig_configs = normalize_sig_configs(sig_configs, 'energy');
catch
    disp('Warning: normalize_sig_configs not found or failed, skipping normalization.');
end


%% 2. Run Simulations
% RLC parameters
omega0 = 1.0;   % Tuned exactly to the local superoscillatory frequency
L = 1;
C0 = 1 / (omega0^2 * L); 
R = 0.25;       % Lower damping for more pronounced localized resonance
alpha = 0.15;   % Slightly higher nonlinearity for a more pronounced 2nd harmonic

tspan = [-500, 500];
y0 = [0; 0]; % initial charge and current

for i = 1:length(sig_configs)
    % The signal from generator is complex, we take the real part to drive the RLC circuit
    V_in = @(t) real(sig_configs(i).data(t));
    
    [tgrid, q, I, V_c] = solve_varactor(L, R, C0, alpha, V_in, tspan, y0);
    
    sig_configs(i).tgrid = tgrid;
    sig_configs(i).V_in = V_in(tgrid);
    sig_configs(i).V_c = V_c;
    sig_configs(i).I = I; % Output current for FFT
end

%% 3. Plot Time-Domain Voltages

% Figure 2: Output Voltages
figure;
ax2 = gca;
hold on;
h_out = gobjects(1, length(sig_configs));
for i = 1:length(sig_configs)
    h_out(i) = plot(sig_configs(i).tgrid, sig_configs(i).V_c, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end
grid on;
xlabel('Time [arb]');
ylabel('Output Voltage $V_{c}(t)$', 'Interpreter', 'latex');
title('Nonlinear Capacitor Output Voltage');
PlotUtils.styleAxes(ax2); 

%% Plot Signals

plot_signals(sig_configs);


%% 4. Plot Windowed Output FFT (Super-Region Analysis)
figure;
hold on;

dt_uniform = 0.02; % High resolution for harmonics
N_fft = 2^16;      % Large N for fine frequency resolution
f_sampling_out = 1 / dt_uniform;
freq_axis_out = linspace(-f_sampling_out/2, f_sampling_out/2, N_fft) * 2 * pi;

% Define the super-region window
t_window = [-5, 16]; 

h_arr_out_fft = gobjects(1, length(sig_configs));
for i = 1:length(sig_configs)
    % Resample the output current uniformly
    t_uniform = min(sig_configs(i).tgrid) : dt_uniform : max(sig_configs(i).tgrid);
    I_uniform = interp1(sig_configs(i).tgrid, sig_configs(i).I, t_uniform, 'spline');
    
    % Isolate the super-region
    idx_win = find(t_uniform >= t_window(1) & t_uniform <= t_window(2));
    t_win = t_uniform(idx_win);
    I_win = I_uniform(idx_win);
    
    % Apply a Hann window to smoothly truncate and isolate the behavior
    I_win_smoothed = I_win .* hann(length(I_win))';
    
    % Compute Windowed FFT magnitude
    I_fft = fftshift(abs(fft(I_win_smoothed, N_fft))) / length(t_win);
    
    % Plot it
    h_arr_out_fft(i) = plot(freq_axis_out, I_fft, '-', 'Color', sig_configs(i).color, 'LineWidth', 1.5, 'DisplayName', sig_configs(i).name);
end
legend(h_arr_out_fft, 'Interpreter', 'latex');

% Zoom in on the relevant angular frequencies
xlim([0, 3]);
xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$', 'Interpreter', 'latex');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$', 'Interpreter', 'latex');
title('Windowed Output Current Spectrum (t \in [-15, 15])');

try
    PlotUtils.styleAxes(gca);
catch
end
hold off;

%% 5. Plot Output Spectrogram (Short-Time Fourier Transform)
% Visualizing the localized high-frequency behavior
for i = 1:length(sig_configs)
    if contains(sig_configs(i).name, 'SO')
        figure;
        
        t_uniform = min(sig_configs(i).tgrid) : dt_uniform : max(sig_configs(i).tgrid);
        I_uniform = interp1(sig_configs(i).tgrid, sig_configs(i).I, t_uniform, 'spline');
        
        fs = 1 / dt_uniform;
        window_length = round(15 / dt_uniform);
        noverlap = round(14.5 / dt_uniform); % High overlap for smooth time resolution
        nfft = 1024;
        
        [s, f, t_spec] = spectrogram(I_uniform, hamming(window_length), noverlap, nfft, fs);
        omega_spec = f * 2 * pi;
        t_spec_shifted = t_spec + min(sig_configs(i).tgrid);
        
        % Focus on omega around 2 and time around super-region
        idx_w = find(omega_spec > 0 & omega_spec < 3);
        idx_t = find(t_spec_shifted > -30 & t_spec_shifted < 30);
        
        imagesc(t_spec_shifted(idx_t), omega_spec(idx_w), 20*log10(abs(s(idx_w, idx_t))));
        axis xy;
        colormap jet;
        c = colorbar;
        ylabel(c, 'Magnitude [dB]');
        xlabel('Time [arb]');
        ylabel('Angular Frequency \omega [\omega_0]');
        title('Spectrogram of SO Output Current (Localized Resonance)');
        try PlotUtils.styleAxes(gca); catch; end
    end
end


%% 6. High-Resolution Autoregressive (Burg) Spectral Estimation
% The standard FFT is fundamentally limited by the uncertainty principle: 
% a short time window inherently causes spectral smearing (broad peaks).
% To resolve sharp peaks in a short window, we use a parametric AR model.
figure;
hold on;

p_order = 30; % Autoregressive model order

h_arr_ar = gobjects(1, length(sig_configs));
for i = 1:length(sig_configs)
    t_uniform = min(sig_configs(i).tgrid) : dt_uniform : max(sig_configs(i).tgrid);
    I_uniform = interp1(sig_configs(i).tgrid, sig_configs(i).I, t_uniform, 'spline');
    idx_win = find(t_uniform >= t_window(1) & t_uniform <= t_window(2));
    I_win = I_uniform(idx_win);
    
    % Compute AR spectrum
    [pxx, w] = pburg(I_win, p_order, N_fft, f_sampling_out);
    w_angular = w * 2 * pi;
    
    % Normalize for plotting
    pxx = pxx / max(pxx(w_angular > 0.1 & w_angular < 3));
    
    h_arr_ar(i) = plot(w_angular, 10*log10(pxx), '-', 'Color', sig_configs(i).color, 'LineWidth', 1.5, 'DisplayName', sig_configs(i).name);
end
legend(h_arr_ar, 'Interpreter', 'latex');

xlim([0, 3]);
ylim([-60, 5]);
xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$', 'Interpreter', 'latex');
ylabel('\boldmath$\mathbf{PSD \ Estimate \ [dB]}$', 'Interpreter', 'latex');
title('High-Resolution AR Spectral Estimation (Burg)');

try PlotUtils.styleAxes(gca); catch; end
hold off;
