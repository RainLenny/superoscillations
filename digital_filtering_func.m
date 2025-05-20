function digital_filtering_func(syms_signal,syms_superoscillation,angular_freqs, syms_filter, f_sampling,window_length_for_filter)
%% SAMPLING
% fundamental_period of the superoscillating signal:
T_period = compute_fundamental_period(angular_freqs,f_sampling); 
signals_duration = T_period * 4; % total duration to simulate
dt = 1 / f_sampling; % sample‐interval in seconds
t_axis = -signals_duration/2 : dt : signals_duration/2;

% Sample the symbolic signal
sample = @(expr) double( subs(expr, symvar(expr), t_axis) );
% now sample each symbolic signal in one line:
[sampled_signal, sampled_superoscillation, sampled_filter] = deal( ...
    sample(syms_signal), ...
    sample(syms_superoscillation), ...
    sample(syms_filter)  ...
);

%% Windowing filter
% Find indices where the condition t ∈ [-window_length_for_filter, window_length_for_filter] is met
valid_indices = abs(t_axis) <= window_length_for_filter;

% Zero out samples outside the range
sampled_filter(~valid_indices) = 0;

% Remove zeroed samples (keep only non-zero samples)
sampled_filter = sampled_filter(valid_indices);
filtered_t_axis = t_axis(valid_indices);

%% Filtering
filtered_signal = conv(sampled_signal, sampled_filter, 'same')*dt;
filtered_superoscilation = conv(sampled_superoscillation, sampled_filter, 'same')*dt;

%% Plots

%% Plot signals and filter
figure
plot(t_axis, real(sampled_superoscillation), '-o','color', 'b', 'LineWidth', 4, 'DisplayName', 'Superoscillation');
hold on;
plot(t_axis, sampled_signal, '-o','color', 'r', 'LineWidth', 4, 'DisplayName', 'Sinc');
plot(filtered_t_axis, sampled_filter, '-o','color', 'k', 'LineWidth', 4, 'DisplayName', 'Filter');
xlabel('Time (s)');
ylabel('Amplitude');
legend('show');
grid on;

%% Plot FFT of signals and filter
% Compute FFT
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

fft_superoscilation = fftshift(abs(fft(sampled_superoscillation,N)));
fft_sinc = fftshift(abs(fft(sampled_signal,N)));
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
xlim([-2,2])

%% Plot filtered signals
figure;

% Create tiled layout with shared axis labels
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Subplot 1: Superoscillating signal
ax1 = nexttile;
plot(t_axis, real(sampled_superoscillation), 'r', 'LineWidth', 6, 'DisplayName', 'Original');
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
plot(t_axis, sampled_signal, 'color', 'r', 'LineWidth', 6, 'DisplayName', 'Original');
hold on;
plot(t_axis, filtered_signal, ':b', 'LineWidth', 6, 'DisplayName', 'Filtered');
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
end

