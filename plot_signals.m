function [] = plot_signals(freq_scaling,amp_scaling)
%% Importing signals
[VinSOFun, VinCosFun, angular_freqs,angular_freqs_COS] = generate_signals_for_plot(freq_scaling,amp_scaling);


%% SAMPLING Constants
f_sampling = 60/(2*pi);
% fundamental_period of the superoscillating signal:
T_period = compute_fundamental_period([angular_freqs,angular_freqs_COS],f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample‐interval in seconds
t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1);
t_axis = t_axis(:);


%% Sample the signals
sampled_signal = VinCosFun(t_axis);
sampled_superoscillation = VinSOFun(t_axis);

%% Plot in time
figure
hold on;
plot(t_axis, real(sampled_signal), '-','color', 'b', 'LineWidth', 4, 'DisplayName', 'COS');
plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'LineWidth', 4, 'DisplayName', 'SO');

xlabel('Time [arb]', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Amplitude [arb]', 'FontWeight', 'bold','FontSize',12);
legend('show', 'FontWeight', 'bold','FontSize',12);
grid off;
xlim([100,400])
ax=gca;
ax.FontWeight = 'bold'
ax.FontSize = 13

%% FFT the signals
% Compute FFT
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N)))*dt;
fft_superoscillation = fft_superoscillation/(sum(fft_superoscillation));

fft_cos = fftshift(abs(fft(sampled_signal,N)))*dt;
fft_cos = fft_cos/(sum(fft_cos));
% fft_cos = fftshift(abs(fft(sampled_signal,N)))*dt;


%% Plot FFTs
figure;
hold on;
plot(freq_axis, (fft_cos), '-', 'Color', 'b', 'LineWidth', 4, 'DisplayName', 'COS');

plot(freq_axis, fft_superoscillation, '-', 'Color', 'r', 'LineWidth', 4, 'DisplayName', 'SO');

% Customize plot
xlabel('Angular frequency [arb]', 'FontWeight', 'bold','FontSize',12);
ylabel('Amplitude [arb]', 'FontWeight', 'bold','FontSize',12);
legend('show', 'FontWeight', 'bold','FontSize',12,'Location', 'Best');
grid off;
xlim([-1,1])
% ylim([0,0.5])
ax=gca;
ax.FontWeight = 'bold'
ax.FontSize = 13
end

