%% Importing signals
[VinSOFun, VinCosFun, angular_freqs,angular_freqs_COS] = generate_signals_for_plot_no_gaussian(1,1);


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
cos1 = -cos(t_axis);
sampled_superoscillation = VinSOFun(t_axis);

%% Plot in time
figure
hold on;
plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'LineWidth', 4, 'DisplayName', '\textbf{SO}');
plot(t_axis, cos1, '--','color', 'black', 'LineWidth', 4, 'DisplayName', '\boldmath$\mathbf{\omega_0}$');

% Updated X-Label with bold math and text
xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$','FontSize', 14, 'Interpreter', 'latex');
legend('show', 'FontWeight', 'bold','FontSize',14,'Location', 'best','Interpreter', 'latex');
grid off;
xlim([0,25])
ax=gca;
ax.FontWeight = 'bold'
ax.FontSize = 13

ax.TickLabelInterpreter = 'latex';

ax.FontWeight = 'bold';

xticks = ax.XTick;
xticklabels = arrayfun(@(x) sprintf('$\\mathbf{%g}$', x), xticks, 'UniformOutput', false);
ax.XTickLabel = xticklabels;

yticks = ax.YTick;
yticklabels = arrayfun(@(y) sprintf('$\\mathbf{%g}$', y), yticks, 'UniformOutput', false);
ax.YTickLabel = yticklabels;

%% Plot in time
figure
hold on;
% Capture the handles (h1, h2) as you plot
h2 = plot(t_axis, real(sampled_signal), '-','color', 'b', 'LineWidth', 4, 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'LineWidth', 4, 'DisplayName', '\textbf{SO}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$','FontSize', 14, 'Interpreter', 'latex');

% Pass the handles in the specific order you want them to appear
legend([h1, h2], 'FontWeight', 'bold', 'FontSize', 14, 'Location', 'best', 'Interpreter', 'latex');

grid off;
xlim([0,50])
ax=gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;

ax.TickLabelInterpreter = 'latex';

ax.FontWeight = 'bold';

xticks = ax.XTick;
xticklabels = arrayfun(@(x) sprintf('$\\mathbf{%g}$', x), xticks, 'UniformOutput', false);
ax.XTickLabel = xticklabels;

yticks = ax.YTick;
yticklabels = arrayfun(@(y) sprintf('$\\mathbf{%g}$', y), yticks, 'UniformOutput', false);
ax.YTickLabel = yticklabels;

% %% FFT the signals
% % Compute FFT
% N = length(t_axis);
% freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis
% 
% fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N)))*dt;
% fft_superoscillation = fft_superoscillation/(sum(fft_superoscillation));
% 
% fft_cos = fftshift(abs(fft(sampled_signal,N)))*dt;
% fft_cos = fft_cos/(sum(fft_cos));
% % fft_cos = fftshift(abs(fft(sampled_signal,N)))*dt;
% 
% 
% %% Plot FFTs
% figure;
% hold on;
% plot(freq_axis, (fft_cos), '-s', 'Color', 'b', 'LineWidth', 4, 'DisplayName', 'COS');
% 
% plot(freq_axis, fft_superoscillation, '-o', 'Color', 'r', 'LineWidth', 4, 'DisplayName', 'SO');
% 
% % Customize plot
% xlabel('Angular frequency [arb]', 'FontWeight', 'bold','FontSize',12);
% ylabel('Amplitude [arb]', 'FontWeight', 'bold','FontSize',12);
% legend('show', 'FontWeight', 'bold','FontSize',12,'Location', 'Best');
% grid off;
% xlim([-1.2,1.2])
% % ylim([0,0.5])
% ax=gca;
% ax.FontWeight = 'bold'
% ax.FontSize = 13
