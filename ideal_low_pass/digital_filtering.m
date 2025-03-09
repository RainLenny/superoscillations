clear all
generate_signals

%% CONSTANTS
f_sampling = 15;
corner_freq_of_filter = 0.7;
% window_length_for_filter = pi/0.7;
window_length_for_filter = inf;

% Compute the GCD of the frequencies
gcd_omega = angular_freqs(1);
for i = 2:length(angular_freqs)
    gcd_omega = gcd(sym(gcd_omega), sym(angular_freqs(i)));  % Use symbolic for precision
end

% Compute the fundamental period
T_period = 2*pi / double(gcd_omega);

%% Sampling
signals_duration = T_period*4; % In seconds

t_axis = -signals_duration/2:1/f_sampling:signals_duration/2; % Time samples
dt = 1/f_sampling;

% Sample the symbolic signal
sampled_superoscilation = double(subs(superoscillations_signal, t, t_axis));
sampled_sinc = double(subs(sinc_signal, t, t_axis));
sampled_filter = double(subs(filter_signal, t, t_axis));

%% Windowing filter
% Find indices where the condition t ∈ [-window_length_for_filter, window_length_for_filter] is met
valid_indices = abs(t_axis) <= window_length_for_filter;

% Zero out samples outside the range
sampled_filter(~valid_indices) = 0;

% Remove zeroed samples (keep only non-zero samples)
sampled_filter = sampled_filter(valid_indices);
filtered_t_axis = t_axis(valid_indices);

%% Filtering
filtered_sinc = conv(sampled_sinc, sampled_filter, 'same')*dt;
filtered_superoscilation = conv(sampled_superoscilation, sampled_filter, 'same')*dt;

%% Plots

%% Plot signals and filter
figure
plot(t_axis, real(sampled_superoscilation), '-o','color', 'b', 'LineWidth', 4, 'DisplayName', 'Superoscillation');
hold on;
plot(t_axis, sampled_sinc, '-o','color', 'r', 'LineWidth', 4, 'DisplayName', 'Sinc');
plot(filtered_t_axis, sampled_filter, '-o','color', 'k', 'LineWidth', 4, 'DisplayName', 'Filter');
xlabel('Time (s)');
ylabel('Amplitude');
legend('show');
grid on;

%% Plot FFT of signals and filter
% Compute FFT
N = length(t_axis)+10000;
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

fft_superoscilation = fftshift(abs(fft(sampled_superoscilation,N)));
fft_sinc = fftshift(abs(fft(sampled_sinc,N)));
fft_filter = fftshift(abs(fft(sampled_filter,N)));

% Plot FFTs
figure;
plot(freq_axis, fft_sinc, '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Sinc Signal');
hold on;
plot(freq_axis, fft_superoscilation, '-o', 'Color', 'g', 'LineWidth', 2, 'DisplayName', 'Superoscillation Signal');
plot(freq_axis, fft_filter, '-o', 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Filter Signal');

% Customize plot
xlabel('Frequency (Rad/s)');
ylabel('Amplitude');
title('FFT of Signals');
legend('show');
grid on;

%% Plot filtered signals
figure;

% Create tiled layout with shared axis labels
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Subplot 1: Superoscillating signal
ax1 = nexttile;
plot(t_axis, real(sampled_superoscilation), 'r', 'LineWidth', 6, 'DisplayName', 'Original');
hold on;
plot(t_axis, real(filtered_superoscilation), ':b', 'LineWidth', 6, 'DisplayName', 'Filtered');
title('Superoscillating Signal');
legend('show', 'fontweight', 'bold', 'fontsize', 12);
grid on;
xlim([-15, 15]);
ax1.FontWeight = 'bold';
ax1.FontSize = 12;
ylim([-0.5, 1]);

% Subplot 2: Sinc signal
ax2 = nexttile;
plot(t_axis, sampled_sinc, 'color', 'r', 'LineWidth', 6, 'DisplayName', 'Original');
hold on;
plot(t_axis, filtered_sinc, ':b', 'LineWidth', 6, 'DisplayName', 'Filtered');
title('Sinc Signal');
legend('show', 'fontweight', 'bold', 'fontsize', 12);
grid on;
xlim([-15, 15]);
ax2.FontWeight = 'bold';
ax2.FontSize = 12;
ylim([-0.5, 1]);

% Shared axis labels
xlabel(t, ['Time [2' char(960) '/' char(969) ']'], 'FontWeight', 'bold', 'FontSize', 12);
ylabel(t, 'Amplitude [arb. units]', 'FontWeight', 'bold', 'FontSize', 12);

hold off;