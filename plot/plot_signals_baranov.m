%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();
[VinSOFun_comp, angular_freqs] = generate_SO_Baranov(1,1,true);
[VinCosFun_comp_unscaled, angular_freqs_COS] = generate_Cos_reference(0.7, -1, true);

% Scale Cos_signal symbolically to match the peak of SO_signal
[~, peak_SO] = normalize_signals(VinSOFun_comp, 'peak');
[~, peak_COS] = normalize_signals(VinCosFun_comp_unscaled, 'peak');
scale_factor = peak_SO / peak_COS;

[VinCosFun_comp, angular_freqs_COS] = generate_Cos_reference(0.7, -scale_factor, true);

% Normalize both signals symbolically by the peak of the first signal
[VinSOFun_comp, VinCosFun_comp] = normalize_signals({VinSOFun_comp, VinCosFun_comp}, 'peak');

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

xlim([100,400]);
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
h2 = plot(freq_axis, fft_cos, '-', 'Color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
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

% 5. Axis Formatting (Scales y by 1e3 and adds top-left exponent \times 10^{-3})
PlotUtils.styleAxes(gca, 1e3);
hold off;

%% Calculate ratio between peak amplitude and amplitude at resonance 
%Stable Ratio Calculation
% 1. Find the peak (Highest Amplitude)
[max_amp, max_idx] = max(fft_superoscillation);
freq_at_max = freq_axis(max_idx);

% 2. USE INTERPOLATION for Frequency 1.0
% This removes the dependency on N by finding the value on the curve at exactly 1.0
target_freq = 1.0;
amp_at_1 = interp1(freq_axis, fft_superoscillation, target_freq, 'pchip');

% 3. Calculate the ratio
amp_ratio = max_amp / amp_at_1;

%Display Results
fprintf('\n--- Stable Amplitude Ratio Analysis ---\n');
fprintf('N (FFT points):            %d\n', N);
fprintf('Peak Frequency Found:      %.8f omega_0\n', freq_at_max);
fprintf('Target Frequency:          %.8f omega_0\n', target_freq);
fprintf('----------------------------------------\n');
fprintf('Highest Amplitude (Max):   %.6e\n', max_amp);
fprintf('Amplitude at exactly 1.0:  %.6e\n', amp_at_1);
fprintf('Ratio (Max / Amp@1):       %.6e\n', amp_ratio);
fprintf('----------------------------------------\n');