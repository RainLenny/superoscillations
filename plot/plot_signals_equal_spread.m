%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();
[VinSOFun_comp, angular_freqs] = generate_SO_equal_spread(1,1,true);
angular_freqs_COS = angular_freqs;
t_0 = 250; T = 100;
s_t_COS_unscaled = @(t) sum(1 .* exp(1i * angular_freqs_COS .* t(:)), 2);
Cos_signal_unscaled_for_norm = @(t) reshape(conj(s_t_COS_unscaled(t) .* exp(-(t(:)-t_0).^2./T^2)), size(t));

% Scale Cos_signal symbolically to match the energy of SO_signal
[~, n_factor_SO] = normalize_signals(VinSOFun_comp, 'energy');
[~, n_factor_COS] = normalize_signals(Cos_signal_unscaled_for_norm, 'energy');
scale_factor = n_factor_SO / n_factor_COS;

[VinCosFun_comp, angular_freqs_COS] = generate_Cos_reference(angular_freqs_COS, scale_factor, true);

% Normalize both signals symbolically by the energy of the first signal
[VinSOFun_comp, VinCosFun_comp] = normalize_signals({VinSOFun_comp, VinCosFun_comp}, 'energy');

VinSOFun = @(t) real(VinSOFun_comp(t));
VinCosFun = @(t) real(VinCosFun_comp(t));


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

figure
hold on;
% Capture the handles (h1, h2) as you plot
h2 = plot(t_axis, real(sampled_signal), '-','color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'DisplayName', '\textbf{SO}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');

% Pass the handles in the specific order you want them to appear
legend([h1, h2]);

PlotUtils.styleAxes(gca);
hold off;

%% FFT the signals
% Compute FFT
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

fft_cos = fftshift(abs(fft(sampled_signal,N))) / N;

fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N))) / N;


%% Plot FFTs
figure;
hold on;

% 1. Plotting Data
h2 = plot(freq_axis, (fft_cos), '-', 'Color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(freq_axis, fft_superoscillation, '-', 'Color', 'r', 'DisplayName', '\textbf{SO}');

% 2. Vertical reference line & Annotation
xline(1.0, 'Color', 'black', 'LineWidth', 4);
text(1, 4.5e-3, '\boldmath$\mathbf{\omega_0}$', ...
    'Color', 'k', 'FontSize', 18, 'Rotation', 90, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% 3. Standard Labels
xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');

% 4. Legend and Limits
legend([h1, h2]);
xlim([0, 1.1]);
ylim([0, 4e-2]);

PlotUtils.styleAxes(gca);
hold off;

