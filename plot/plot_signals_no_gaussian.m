clear; clc;
%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();
[VinSOFun_comp, angular_freqs] = generate_SO_from_dat_file('SO_Baranov', 1, 1, false);
[VinCosFun_comp, angular_freqs_COS] = generate_Cos_reference(0.9, -1, false);



% Normalize both signals symbolically by the peak of the first signal (using version with Gaussian envelope)
[VinSOFun_comp,VinCosFun_comp] = normalize_signals({VinSOFun_comp,VinCosFun_comp}, 'energy');



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
cos1 = -cos(t_axis);
sampled_superoscillation = VinSOFun(t_axis);

%% Plot in time
figure
hold on;
plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'DisplayName', '\textbf{SO}');
plot(t_axis, cos1,'color', 'black', 'DisplayName', '\boldmath$\mathbf{\omega_0}$');

% Updated X-Label with bold math and text
xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');
legend('show');
xlim([-15,30]);

PlotUtils.styleAxes(gca);
hold off;


%% Plot in time
figure
hold on;
% Capture the handles (h1, h2) as you plot
h2 = plot(t_axis, real(sampled_signal), '-','color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'DisplayName', '\textbf{SO}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');

% Pass the handles in the specific order you want them to appear
legend([h1, h2]);
xlim([-15,30]);

PlotUtils.styleAxes(gca);
hold off;


%% FFT the signals
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

% Calculate FFT and normalize by N for true amplitude
fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N))) / N;

fft_cos = fftshift(abs(fft(sampled_signal,N))) / N;




%% Plot FFTs
figure;
hold on;
% Updated DisplayNames to match the bold styling
h1 = plot(freq_axis, fft_cos, '-o', 'Color', 'b','DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');

h2 = plot(freq_axis, fft_superoscillation, '-o', 'Color', 'r', 'DisplayName', '\textbf{SO}');

% 2. Vertical reference line & Annotation
xline(1.0, 'Color', 'black','LineWidth', 6);
text(1, 30e-2, '\boldmath$\mathbf{\omega_0}$', ...
    'Color', 'k', 'FontSize', 18, 'Rotation', 90, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% Updated X and Y labels with bold math/text formatting
xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');

legend([h2,h1]);
grid off;
xlim([0,1.1])
% ylim([0,0.5])

PlotUtils.styleAxes(gca);
hold off;