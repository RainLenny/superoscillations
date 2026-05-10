%% 1. Initialization and Parameters
clear; clc; close all;

freq_scaling = 1;
amp_scaling = 1;

% SO Parameters
angular_freqs_SO = [1,2,3,4,5] * 0.18 * freq_scaling;
amps_SO    = [-0.156067704462866 + 0.331660754319902i,...
    -0.861836830340772 - 1.041781767817215i,...
    2.340666531884434 - 0.600981019561177i,...
    -0.502399901009243 + 2.633672512223463i,...
    -1.820362096071554 - 1.322570479164972i] * amp_scaling;

% Reference Parameters
angular_freqs_COS = angular_freqs_SO; 
Cos_normalization = 2;

% Function handles
% (Corrected Cos_signal to use angular_freqs_COS)
Cos_signal = @(t)  conj(sum( Cos_normalization .* exp(1i * angular_freqs_COS .* t),2));
SO_signal = @(t)  conj(sum( amps_SO .* exp(1i * angular_freqs_SO .* t),2));

%% 2. Sampling Constants & Time Axis Definition
f_sampling = 60/(2*pi);

% Assume 'compute_fundamental_period' is defined elsewhere in your path
T_period = compute_fundamental_period([angular_freqs_SO, angular_freqs_COS], f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample-interval in seconds

t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1);
t_axis = t_axis(:);

% Sample the global signals
sampled_signal = Cos_signal(t_axis);
sampled_superoscillation = SO_signal(t_axis);


%% Plot in time
figure
hold on;
% Capture the handles (h1, h2) as you plot
h2 = plot(t_axis, real(sampled_signal), '-','color', 'b', 'LineWidth', 4, 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'LineWidth', 4, 'DisplayName', '\textbf{SO}');

xlim([-50,50])
%% 3. Truncate Signals at Independent Time Windows
% A. Define the window for the Superoscillating (SO) ripples
t_start_so = 0; 
t_end_so = 4*pi;
idx_so = (t_axis >= t_start_so) & (t_axis <= t_end_so);
so_trunc = sampled_superoscillation(idx_so);
t_trunc_so = t_axis(idx_so);

% B. Define the window for the Reference ripples
% Adjust t_start_ref to wherever the blue line has nice, uniform ripples
t_start_ref = 10.9956; 
t_end_ref = t_start_ref + (t_end_so - t_start_so); % Force the same duration
idx_ref = (t_axis >= t_start_ref) & (t_axis <= t_end_ref);
cos_trunc = sampled_signal(idx_ref);
t_trunc_ref = t_axis(idx_ref);

% C. Enforce exactly equal array lengths (to prevent dt rounding mismatches)
min_len = min(length(so_trunc), length(cos_trunc));
so_trunc = so_trunc(1:min_len);
cos_trunc = cos_trunc(1:min_len);
t_trunc_so = t_trunc_so(1:min_len);
t_trunc_ref = t_trunc_ref(1:min_len);

% D. Apply a Window Function (Hamming) to prevent spectral leakage
window = hamming(min_len);
so_trunc_win = so_trunc .* window;
cos_trunc_win = cos_trunc .* window;

%% 4. Plot the Truncated Time-Domain Signals 
% Shift the time axes to 0 just for this plot so we can compare wave shapes directly
t_shifted = t_trunc_so - t_trunc_so(round(min_len/2)); 

figure;
hold on;
plot(t_shifted, real(cos_trunc)+0.9549032, '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', '\boldmath$\mathbf{0.9\omega_0 \ (Ref \ Window)}$');
plot(t_shifted, real(so_trunc), '-', 'Color', 'r', 'LineWidth', 3, 'DisplayName', '\textbf{SO \ (Ripple \ Window)}');
xlabel('\boldmath$\mathbf{Relative \ Time \ [2\pi/\omega_0]}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$', 'FontSize', 14, 'Interpreter', 'latex');
title('\textbf{Isolated Ripple Regions (Time-Shifted for Comparison)}', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontWeight', 'bold', 'FontSize', 14, 'Location', 'best', 'Interpreter', 'latex');
grid on;

%% 5. Compute Local FFT
% Because both slices are the exact same length, we use one frequency axis
freq_axis_trunc = linspace(-f_sampling/2, f_sampling/2, min_len) * 2 * pi;

% Calculate FFT (using windowed signals)
fft_so_trunc = fftshift(abs(fft(so_trunc_win, min_len))) / min_len;
fft_cos_trunc = fftshift(abs(fft(cos_trunc_win, min_len))) / min_len;

%% 6. Plot the Local Spectrum
figure;
hold on;
plot(freq_axis_trunc, fft_cos_trunc, '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
plot(freq_axis_trunc, fft_so_trunc, '-', 'Color', 'r', 'LineWidth', 3, 'DisplayName', '\textbf{SO}');

xlim([-3, 3]); % Zoom in on the relevant frequency range to see the shift
xlabel('\boldmath$\mathbf{Local \ Angular \ frequency \ [\omega_0]}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('\boldmath$\mathbf{Magnitude}$', 'FontSize', 14, 'Interpreter', 'latex');
title('\textbf{Local Spectrum Comparison (FFT of Truncated Windows)}', 'FontSize', 14, 'Interpreter', 'latex');
legend('FontWeight', 'bold', 'FontSize', 14, 'Location', 'best', 'Interpreter', 'latex');
grid on;