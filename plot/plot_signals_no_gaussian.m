%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();
[VinSOFun_comp, angular_freqs] = generate_SO_from_dat_file('SO_Baranov', 1, 1, false);
[VinSOFun_comp_for_norm, ~] = generate_SO_from_dat_file('SO_Baranov', 1, 1, true);
[VinCosFun_comp_unscaled, angular_freqs_COS] = generate_Cos_reference(0.7, -1, false);
[VinCosFun_comp_for_norm, ~] = generate_Cos_reference(0.7, -1, true);

% Scale Cos_signal symbolically to match the peak of SO_signal
[~, peak_SO] = normalize_signals(VinSOFun_comp_for_norm, 'peak');
[~, peak_COS] = normalize_signals(VinCosFun_comp_for_norm, 'peak');
scale_factor = peak_SO / peak_COS;

[VinCosFun_comp, angular_freqs_COS] = generate_Cos_reference(0.7, -scale_factor, false);

% Normalize both signals symbolically by the peak of the first signal (using version with Gaussian envelope)
[~, N_factor] = normalize_signals(VinSOFun_comp_for_norm, 'peak');

VinSOFun_comp = @(t) VinSOFun_comp(t) / N_factor;
VinCosFun_comp = @(t) VinCosFun_comp(t) / N_factor;

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
plot(t_axis, cos1, '--','color', 'black', 'DisplayName', '\boldmath$\mathbf{\omega_0}$');

% Updated X-Label with bold math and text
xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');
legend('show');
xlim([0,25]);

PlotUtils.styleAxes(gca);
hold off;

%% Plot in time
figure
hold on;
% Capture the handles (h1, h2) as you plot
h2 = plot(t_axis, real(sampled_signal), '-','color', 'b', 'DisplayName', '\boldmath$\mathbf{0.9\omega_0}$');
h1 = plot(t_axis, real(sampled_superoscillation), '-','color', 'r', 'DisplayName', '\textbf{SO}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');

% Pass the handles in the specific order you want them to appear
legend([h1, h2]);
xlim([0,50]);

PlotUtils.styleAxes(gca);
hold off;
