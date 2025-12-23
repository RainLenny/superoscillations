function [] = plot_signals(freq_scaling,amp_scaling)
%% Importing signals
[VinSOFun, VinCosFun, angular_freqs,angular_freqs_COS] = generate_signals(freq_scaling,amp_scaling);


%% SAMPLING Constants
f_sampling = 60/(2*pi);
% fundamental_period of the superoscillating signal:
T_period = compute_fundamental_period([angular_freqs,angular_freqs_COS],f_sampling);
signals_duration = T_period * 10; % total duration to simulate
dt = 1 / f_sampling; % sample‐interval in seconds
t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1);
t_axis = t_axis(:);


%% Sample the signals
sampled_signal = VinCosFun(t_axis);
sampled_superoscillation = VinSOFun(t_axis);

%% Plot in time
figure
plot(t_axis, real(sampled_superoscillation), '-o','color', 'b', 'LineWidth', 4, 'DisplayName', 'Superoscillation');
hold on;
plot(t_axis, real(sampled_signal), '-o','color', 'r', 'LineWidth', 4, 'DisplayName', 'Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('show');
grid on;
% xlim([-25,25])

%% FFT the signals
% Compute FFT
N = length(t_axis)+1000;
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N)))*dt;
fft_superoscillation = fft_superoscillation/(max(fft_superoscillation));
fft_cos = fftshift(abs(fft(sampled_signal,N)))*dt;


%% Plot FFTs
figure;
plot(freq_axis, fft_cos/1800, '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Signal');
hold on;
plot(freq_axis, fft_superoscillation, '-o', 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Superoscillations');
% Customize plot
xlabel(['Frequency [' char(969) ']'], 'FontWeight', 'bold','FontSize',12);
ylabel('Amplitude [arb. units]', 'FontWeight', 'bold','FontSize',12);
title('Frequncy response of signals', 'FontWeight', 'bold','FontSize',12);
legend('show', 'FontWeight', 'bold','FontSize',12);
grid on;
xlim([-1.5,1.5])
end

