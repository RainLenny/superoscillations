clear all
generate_signals

%% CONSTANTS
f_sampling = 60/(2*pi);

% Filter Parameters
omega0 = 1*pi+1.5;  % Corner frequency (rad/s)
filter_order = 7;         % Filter order
epsilon = 0.04475;  % Ripple factor 

T_period = compute_fundamental_period(angular_freqs(2:4),f_sampling);


%% Sampling Signal
signals_duration = T_period*20; % In seconds

dt = 1/f_sampling;
t_axis = -signals_duration/2:dt:signals_duration/2; % Time samples
t_axis = t_axis(1:end-1);

% Sample the symbolic signal
sampled_superoscilation = double(subs(superoscillations_signal, t, t_axis));
sampled_cos = double(subs(cos_signal, t, t_axis));

%% Generate Sampled Filter
H_sym = chebychev_filter_magnitude(omega0, filter_order, epsilon);
H_func = matlabFunction(H_sym, 'Vars', sym('omega'));

% angular frequency vector 
domega = 2*pi / signals_duration;
Nfft = length(t_axis);  
% Define angular frequency axis (rad/s):
omega_axis = (-Nfft/2 : Nfft/2 - 1) * (2*pi*f_sampling/Nfft);

% Evaluate H(j*omega) on this frequency grid
H_samples = H_func(omega_axis);
H_samples(find(isnan(H_samples))) = 1;
%inverse FFT
sampled_filter = ifft(fftshift(H_samples));
sampled_filter = fftshift(sampled_filter);
sampled_filter = sampled_filter / sum(sampled_filter);

%% Filtering
filtered_cos = conv(sampled_cos, sampled_filter, 'same');
filtered_superoscilation = conv(sampled_superoscilation, sampled_filter, 'same');

%% Plot signals and filter
figure
plot(t_axis, real(sampled_superoscilation), '-o','color', 'b', 'LineWidth', 4, 'DisplayName', 'Superoscillation');
hold on;
plot(t_axis, sampled_cos, '-o','color', 'r', 'LineWidth', 4, 'DisplayName', 'cos');
plot(t_axis, real(sampled_filter), '-o','color', 'k', 'LineWidth', 4, 'DisplayName', 'Filter');
xlabel('Time (s)');
ylabel('Amplitude');
legend('show');
grid on;

%% Plot fitler before ifft
figure;

% Plot Magnitude Response in dB
subplot(2,1,1);
plot(omega_axis, abs(H_samples), '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Magnitude Response');
hold on;
legend('show');
title('Magnitude Response of the Chebyshev Filter (dB)');
xlabel('Frequency (Rad/s)');
ylabel('Magnitude');
grid on;

% Plot Phase Response
subplot(2,1,2);
plot(omega_axis, unwrap(angle(H_samples)), '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Phase Response');
hold on;
legend('show');
title('Phase Response of the Chebyshev Filter');
xlabel('Frequency (Rad/s)');
ylabel('Phase');
grid on;

%% Plot FFT of signals and filter
% Compute FFT
Nfft = length(t_axis);
omega_axis = linspace(-f_sampling/2, f_sampling/2, Nfft)*2*pi; % Frequency axis

fft_superoscilation = fftshift(abs(fft(sampled_superoscilation,Nfft)));
fft_cos = fftshift(abs(fft(sampled_cos,Nfft)));
fft_filter = fftshift(abs(fft(sampled_filter,Nfft)));

% Plot FFTs
figure;
plot(omega_axis, fft_cos, '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'cos Signal');
hold on;
plot(omega_axis, fft_superoscilation, '-o', 'Color', 'g', 'LineWidth', 2, 'DisplayName', 'Superoscillation Signal');
plot(omega_axis, fft_filter, '-o', 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Filter Signal');

% Customize plot
xlabel('Frequency (Rad/s)','FontWeight', 'bold', 'FontSize', 12);
ylabel('Amplitude [arb. units]','FontWeight', 'bold', 'FontSize', 12);
title('FFT of Signals');
legend('show');
xlim([-8,8])
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
xlim([-6, 6]);
ax1.FontWeight = 'bold';
ax1.FontSize = 12;
ylim([-3, 4]);

% Subplot 2: cos signal
ax2 = nexttile;
plot(t_axis, sampled_cos, 'color', 'r', 'LineWidth', 6, 'DisplayName', 'Original');
hold on;
plot(t_axis, real(filtered_cos), ':b', 'LineWidth', 6, 'DisplayName', 'Filtered');
title('cos Signal');
legend('show', 'fontweight', 'bold', 'fontsize', 12);
grid on;
xlim([-6, 6]);
ax2.FontWeight = 'bold';
ax2.FontSize = 12;
ylim([-1.2, 1.2]);

% Shared axis labels
xlabel(t, ['Time [2' char(960) '/' char(969) ']'], 'FontWeight', 'bold', 'FontSize', 12);
ylabel(t, 'Amplitude [arb. units]', 'FontWeight', 'bold', 'FontSize', 12);

hold off;