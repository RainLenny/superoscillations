% interactive_Eli_signal.m
% Creates an interactive figure to plot the superoscillating signal and its
% instantaneous frequency, as well as its frequency spectrum, with sliders for A, Omega, and delta.

clear; clc; close all;

%% Setup Path
addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
PlotUtils.setupDefaults();

%% Initial Parameters
A_init = 10;
Omega_init = 1.0;
delta_init = 0.1;

% Time axis setup (using constants similar to plot_2p5_cos.m)
f_sampling = 600/(2*pi);
dt = 1 / f_sampling;
t_axis = -100 : dt : 100; 
t_axis = t_axis(:);

%% Create Figure and UI Layout
fig = figure('Name', 'Interactive Eli Signal', 'Position', [100, 100, 900, 800]);

% Create a panel for the plots to leave space at the bottom for sliders
pnl = uipanel('Parent', fig, 'Position', [0.05, 0.25, 0.9, 0.73], 'BorderType', 'none');
tlo = tiledlayout(pnl, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Initial signal computation (no Gaussian envelope for pure freqs)
[signal_func, ~, ~] = generate_signal_Eli(A_init, Omega_init, delta_init, false, false, false);
sig_val_0 = signal_func(0);
if abs(sig_val_0) > 1e-10
    sampled_signal = signal_func(t_axis) / sig_val_0;
else
    sampled_signal = signal_func(t_axis);
end

% Analytical instantaneous frequency calculation
% Formula: Omega + [(A-1)*delta*(A*cos(delta*t) - (A-1))] / [A^2 + (A-1)^2 - 2A(A-1)cos(delta*t)]
num = (A_init - 1) * delta_init .* (A_init .* cos(delta_init .* t_axis) - (A_init - 1));
den = A_init^2 + (A_init - 1)^2 - 2 * A_init * (A_init - 1) .* cos(delta_init .* t_axis);
inst_freq = Omega_init + num ./ den;

% Initial FFT computation (using a much longer time window and Hann window to produce sharp delta-like peaks)
t_fft = (-1000 : dt : 1000).';
if abs(sig_val_0) > 1e-10
    signal_for_fft = signal_func(t_fft) / sig_val_0;
else
    signal_for_fft = signal_func(t_fft);
end
win = hann(length(t_fft));
N_fft = length(t_fft) * 4; % Zero pad
freq_axis = linspace(-f_sampling/2, f_sampling/2, N_fft)*2*pi;
fft_signal = fftshift(abs(fft(signal_for_fft .* win, N_fft)))*dt;
if max(fft_signal) > 1e-10
    fft_signal = 0.9 * fft_signal / max(fft_signal);
end

%% Plotting - Subplot 1 (Time Domain & Inst Freq)
ax1 = nexttile(tlo);
yyaxis(ax1, 'left');
p1 = plot(ax1, t_axis, real(sampled_signal), '-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', '\textbf{SO}');
ylabel(ax1, '\boldmath$\mathbf{Amplitude \ [arb.]}$', 'Interpreter', 'latex');
ax1.YAxis(1).Color = 'k';

yyaxis(ax1, 'right');
p2 = plot(ax1, t_axis, inst_freq, '-', 'Color', '#006400', 'LineWidth', 1.5, 'DisplayName', '\textbf{Inst. Freq.}');
ylabel(ax1, '\boldmath$\mathbf{Inst. \ Freq. \ [\omega_0]}$', 'Interpreter', 'latex');
ax1.YAxis(2).Color = 'k';
% Green text for theoretical max freq (highest Fourier component)
freq_line = yline(ax1, Omega_init, '--', 'Color', '#006400', 'HandleVisibility', 'off', 'LineWidth', 1.5);
freq_text = text(ax1, -23, Omega_init + 0.05, '\boldmath$\mathbf{\Omega}$', 'Color', '#006400', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 14, 'Interpreter', 'latex');

xlabel(ax1, '\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$', 'Interpreter', 'latex');
xlim(ax1, [-27, 27]);
legend(ax1, 'show', 'Location', 'best', 'Interpreter', 'latex');

mask_so = abs(t_axis) <= 5;
update_ylim(ax1, sampled_signal, inst_freq, mask_so, Omega_init);

% Only apply styleAxes initially for baseline styling, we will handle ticks manually in updates
if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(ax1);
end

%% Plotting - Subplot 2 (Frequency Domain)
ax2 = nexttile(tlo);
hold(ax2, 'on');
p3 = plot(ax2, freq_axis, fft_signal, '-', 'Color', 'r', 'LineWidth', 2.5, 'DisplayName', '\textbf{SO}');
xlabel(ax2, '\boldmath$\mathbf{Angular \ Frequency \ [\omega_0]}$', 'Interpreter', 'latex');
ylabel(ax2, '\boldmath$\mathbf{Amplitude \ [scaled]}$', 'Interpreter', 'latex');
xlim(ax2, [0, 1.5]);
ylim(ax2, [0, 1.2]);
legend(ax2, 'show', 'Location', 'best', 'Interpreter', 'latex');

if exist('PlotUtils', 'class')
    PlotUtils.styleAxes(ax2);
end
hold(ax2, 'off');

%% GUI Elements (Sliders)
% Slider A
uicontrol('Parent', fig, 'Style', 'text', 'Position', [50, 100, 120, 20], 'String', 'Amplitude (A)', 'FontSize', 12, 'HorizontalAlignment', 'left');
txt_A = uicontrol('Parent', fig, 'Style', 'text', 'Position', [170, 100, 50, 20], 'String', num2str(A_init), 'FontSize', 12, 'HorizontalAlignment', 'left');
sld_A = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [220, 100, 600, 20], ...
    'Min', 1, 'Max', 100, 'Value', A_init);

% Slider Omega
uicontrol('Parent', fig, 'Style', 'text', 'Position', [50, 65, 120, 20], 'String', 'Omega (\Omega)', 'FontSize', 12, 'HorizontalAlignment', 'left');
txt_Omega = uicontrol('Parent', fig, 'Style', 'text', 'Position', [170, 65, 50, 20], 'String', num2str(Omega_init), 'FontSize', 12, 'HorizontalAlignment', 'left');
sld_Omega = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [220, 65, 600, 20], ...
    'Min', 0.1, 'Max', 5.0, 'Value', Omega_init);

% Slider delta
uicontrol('Parent', fig, 'Style', 'text', 'Position', [50, 30, 120, 20], 'String', 'Delta (\delta)', 'FontSize', 12, 'HorizontalAlignment', 'left');
txt_delta = uicontrol('Parent', fig, 'Style', 'text', 'Position', [170, 30, 50, 20], 'String', num2str(delta_init), 'FontSize', 12, 'HorizontalAlignment', 'left');
sld_delta = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [220, 30, 600, 20], ...
    'Min', 0.001, 'Max', 1.0, 'Value', delta_init);

%% Setup Callbacks
cb = @(es, ed) update_plot(sld_A, sld_Omega, sld_delta, txt_A, txt_Omega, txt_delta, p1, p2, p3, freq_line, freq_text, ax1, ax2, t_axis, t_fft, win, freq_axis, dt, mask_so, N_fft);
addlistener(sld_A, 'ContinuousValueChange', cb);
addlistener(sld_Omega, 'ContinuousValueChange', cb);
addlistener(sld_delta, 'ContinuousValueChange', cb);

%% Helper Functions

function update_plot(sld_A, sld_Omega, sld_delta, txt_A, txt_Omega, txt_delta, p1, p2, p3, freq_line, freq_text, ax1, ax2, t_axis, t_fft, win, freq_axis, dt, mask_so, N_fft)
    % Read slider values
    A = sld_A.Value;
    Omega = sld_Omega.Value;
    delta = sld_delta.Value;
    
    % Update text labels
    txt_A.String = num2str(A, '%.2f');
    txt_Omega.String = num2str(Omega, '%.2f');
    txt_delta.String = num2str(delta, '%.3f');
    
    % Recompute signal (no Gaussian envelope)
    [signal_func, ~, ~] = generate_signal_Eli(A, Omega, delta, false, false, false);
    
    sig_val_0 = signal_func(0);
    if abs(sig_val_0) > 1e-10
        sampled_signal = signal_func(t_axis) / sig_val_0;
    else
        sampled_signal = signal_func(t_axis);
    end
    
    % Analytical instantaneous frequency calculation
    num = (A - 1) * delta .* (A .* cos(delta .* t_axis) - (A - 1));
    den = A^2 + (A - 1)^2 - 2 * A * (A - 1) .* cos(delta .* t_axis);
    inst_freq = Omega + num ./ den;
    
    % Recompute FFT (with zero padding and windowing on the long signal)
    if abs(sig_val_0) > 1e-10
        signal_for_fft = signal_func(t_fft) / sig_val_0;
    else
        signal_for_fft = signal_func(t_fft);
    end
    fft_signal = fftshift(abs(fft(signal_for_fft .* win, N_fft)))*dt;
    if max(fft_signal) > 1e-10
        fft_signal = 0.9 * fft_signal / max(fft_signal);
    end
    
    % Update plots
    p1.YData = real(sampled_signal);
    p2.YData = inst_freq;
    p3.YData = fft_signal;
    
    % Update the theoretical max freq line (highest Fourier component)
    freq_line.Value = Omega;
    freq_text.Position(2) = Omega + 0.05;
    
    % Adjust frequency domain xlim to make sure peaks are visible
    max_freq_to_show = max(1.5, Omega + 0.5);
    xlim(ax2, [0, max_freq_to_show]);
    
    update_ylim(ax1, sampled_signal, inst_freq, mask_so, Omega);
end

function update_ylim(ax, sampled_signal, inst_freq, mask_so, Omega)
    % Reset tick modes to auto to prevent overlapping ticks set by PlotUtils
    ax.YAxis(1).TickValuesMode = 'auto';
    ax.YAxis(1).TickLabelsMode = 'auto';
    ax.YAxis(2).TickValuesMode = 'auto';
    ax.YAxis(2).TickLabelsMode = 'auto';

    if any(mask_so)
        % Dynamic scaling for amplitude
        max_amp = max(abs(real(sampled_signal(mask_so))));
        if max_amp > 0
            yyaxis(ax, 'left');
            ylim(ax, [-max_amp*1.5, max_amp*1.5]);
        end
        
        % Dynamic scaling for frequency (capped to prevent extreme zooming)
        max_freq = max(abs(inst_freq(mask_so)));
        % Cap max_freq relative to Omega so we don't squish everything if there's a singularity
        cap_val = max(5 * Omega, 2.0);
        if max_freq > cap_val
            max_freq = cap_val;
        end
        if max_freq > 0
            yyaxis(ax, 'right');
            ylim(ax, [-max_freq*0.5, max_freq*1.5]);
        end
    end
    
    % Align the 0 of both y-axes by minimally expanding limits
    yyaxis(ax, 'left'); yL = ylim(ax);
    yyaxis(ax, 'right'); yR = ylim(ax);
    if yL(1) < 0 && yL(2) > 0 && yR(1) < 0 && yR(2) > 0
        rL = yL(2) / abs(yL(1));
        rR = yR(2) / abs(yR(1));
        Rs = [rL, rR, 1];
        best_R = Rs(1);
        min_penalty = inf;
        for R = Rs
            new_yL1 = -max(abs(yL(1)), yL(2)/R);
            new_yL2 = max(yL(2), abs(yL(1))*R);
            new_yR1 = -max(abs(yR(1)), yR(2)/R);
            new_yR2 = max(yR(2), abs(yR(1))*R);
            penalty = (new_yL2 - new_yL1)/(yL(2) - yL(1)) + (new_yR2 - new_yR1)/(yR(2) - yR(1));
            if penalty < min_penalty
                min_penalty = penalty;
                best_R = R;
            end
        end
        yyaxis(ax, 'left');  ylim(ax, [-max(abs(yL(1)), yL(2)/best_R), max(yL(2), abs(yL(1))*best_R)]);
        yyaxis(ax, 'right'); ylim(ax, [-max(abs(yR(1)), yR(2)/best_R), max(yR(2), abs(yR(1))*best_R)]);
    end
end
