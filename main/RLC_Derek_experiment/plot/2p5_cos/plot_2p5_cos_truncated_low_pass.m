clear; clc;
%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();
[VinSOFun_comp, angular_freqs, amps_SO] = generate_SO_from_dat_file('2p5_cos_Derek', 1, 1, false);
[VinCosFun_comp, angular_freqs_COS] = generate_Cos_reference(1, 1, false);

VinSOFun = @(t) real(VinSOFun_comp(t)) ./ real(VinSOFun_comp(0));

VinCosFun = @(t) real(VinCosFun_comp(t));


%% SAMPLING Constants
f_sampling = 60/(2*pi);
% fundamental_period of the superoscillating signal:
T_period = compute_fundamental_period([angular_freqs,angular_freqs_COS],f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample-interval in seconds
t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1);
t_axis = t_axis(:);


%% Sample the signals
% Using a cos reference of frequency 1 (\omega_0) as requested
cos1 = VinCosFun(t_axis);
sampled_superoscillation = VinSOFun(t_axis);


%% FFT the signals
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

% Calculate the FFT magnitude
fft_superoscillation = fftshift(abs(fft(sampled_superoscillation, N)))*dt;
fft_cos = fftshift(abs(fft(cos1, N)))*dt;

% % Normalize by the number of samples (N) to get true Fourier coefficients
% fft_superoscillation = fft_superoscillation / (sum(fft_superoscillation));
% fft_cos = fft_cos / (sum(fft_cos));


fft_superoscillation = 0.9 * fft_superoscillation / max(fft_superoscillation);
fft_cos = 0.9 * fft_cos / max(fft_cos);


%% Truncated Filter Definitions
omega_c = 0.7; % Cutoff frequency
T_trunc = pi / 0.7;
% M_func is the multiplier function M(omega) for a frequency omega
M_func = @(omega) (1/pi) * (sinint((omega_c + omega) .* T_trunc) + sinint((omega_c - omega) .* T_trunc));


%% Figure (a): The signals and the truncated filter in time and frequency domains.
figure;
tlo = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Time domain
ax1 = nexttile;
hold on;
% Scale filter to have similar amplitude to signals for visualization
filter_time = (omega_c/pi) * sinc(omega_c * t_axis / pi);
filter_time(abs(t_axis) > T_trunc) = 0; % Truncate in time
filter_time = filter_time / max(filter_time); % scale factor to match image roughly

h_so_time = plot(t_axis, real(sampled_superoscillation), '-', 'color', 'r', 'DisplayName', '\textbf{SO}');
h_cos_time = plot(t_axis, cos1, ':', 'color', 'b', 'DisplayName', '\textbf{\boldmath$\mathbf{\omega_0}$}');
h_filter_time = plot(t_axis, filter_time, '-.', 'color', 'k', 'DisplayName', '\textbf{Filter}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');
xlim([-8, 8]);
PlotUtils.styleAxes(gca);
hold off;


% Subplot 2: Frequency domain
ax2 = nexttile;
hold on;
% Filter in frequency domain: M_func(omega)
filter_freq = M_func(freq_axis);

% Assign handles to each plot
h1 = plot(freq_axis, fft_cos, '-o', 'Color', 'b', 'DisplayName', '\textbf{\boldmath$\mathbf{\omega_0}$}');
h2 = plot(freq_axis, fft_superoscillation, '-o', 'Color', 'r', 'DisplayName', '\textbf{SO}');
h3 = plot(freq_axis, filter_freq , 'Color', 'k', 'DisplayName', '\textbf{Filter}');


xlabel('\boldmath$\mathbf{Frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [scaled]}$');

xlim([0, 1.5]);
ylim([0, 1.2]);
PlotUtils.styleAxes(gca);
hold off;

% Create a single combined legend assigned to the tiledlayout
lgd_all = legend(ax1, [h_filter_time, h_so_time, h_cos_time], '\textbf{Filter}', '\textbf{SO}', '\textbf{\boldmath$\mathbf{\omega_0}$}', 'Orientation', 'horizontal', 'Box', 'off');
lgd_all.Layout.Tile = 'south';


%% Figure (c): The signals filtered by the truncated filter.
% Analytical filtering:
% SO is analytically filtered using the formula on each frequency component
filtered_SO_comp = @(t) real(reshape(conj(sum((amps_SO .* M_func(angular_freqs)) .* exp(1i * angular_freqs .* t(:)), 2)), size(t)));
scale_factor = real(VinSOFun_comp(0));
filtered_superoscillation = filtered_SO_comp(t_axis) ./ scale_factor;

% Cos is analytically filtered (freq 1)
filtered_cos = M_func(angular_freqs_COS) .* cos1;

figure;
% Subplot 1: Superoscillating Signal
subplot(2,1,1);
hold on;
plot(t_axis, real(sampled_superoscillation), '-', 'color', 'r', 'DisplayName', '\textbf{Original}');
plot(t_axis, real(filtered_superoscillation), '--', 'color', 'b', 'DisplayName', '\textbf{Filtered}');
title('\textbf{SO Signal}');
lgd = legend('show', 'Orientation', 'horizontal');
xlim([-12, 12]);
PlotUtils.styleAxes(gca);
% Position legend between the plots without resizing axes
lgd.Units = 'normalized';
drawnow;
lgd.Position(1) = 0.905 - lgd.Position(3);
lgd.Position(2) = 0.46;
hold off;

% Subplot 2: Cos Signal
subplot(2,1,2);
hold on;
plot(t_axis, cos1, '-', 'color', 'r', 'DisplayName', '\textbf{Original}');
plot(t_axis, filtered_cos, '--', 'color', 'b', 'DisplayName', '\textbf{Filtered}');
xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
title('\textbf{\boldmath$\mathbf{COS(\omega_0)}$}');
xlim([-12, 12]);
ylim([-1.2, 1.2]); % Adjust ylim for visibility of the zero line
PlotUtils.styleAxes(gca);
hold off;

% Add shared y-axis label
han = axes('visible', 'off'); 
han.YLabel.Visible = 'on';
ylabel(han, '\boldmath$\mathbf{Amplitude \ [arb.]}$');
PlotUtils.styleAxes(han);
han.Visible = 'off';
han.YLabel.Visible = 'on';


%% Figure (d): The filtered signals on top of each other.
figure;
hold on;

plot(t_axis, real(filtered_superoscillation), '-', 'color', 'r', 'DisplayName', '\textbf{Filtered SO}');
plot(t_axis, filtered_cos, ':', 'color', 'b', 'DisplayName', '\textbf{Filtered \boldmath$\mathbf{\omega_0}$}');

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');
legend('show');
xlim([-6, 6]);
PlotUtils.styleAxes(gca);
hold off;


%% Figure (e): Filtered SO in time with instantaneous frequency
% --- Plot Y-Limit Controls ---
% Set to [min, max] for manual limits, or [] for dynamic scaling
ylim_plot1_so = [-4,7]*10^(-2)/0.0148;
ylim_plot1_freq = [-0.5,1.5];
% -------------------------------

figure;
[inst_freq_num, ~, ~] = compute_instantaneous_frequency(filtered_superoscillation, t_axis);

yyaxis left;
p1 = plot(t_axis, real(filtered_superoscillation), '-', 'Color', 'r', 'DisplayName', '\textbf{Filtered SO}');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');
ax = gca;
ax.YAxis(1).Color = 'k';

% Dynamic scaling for amplitude
mask_so = abs(t_axis) <= 5;
if ~isempty(ylim_plot1_so)
    ylim(ylim_plot1_so);
elseif any(mask_so)
    max_amp = max(abs(real(filtered_superoscillation(mask_so))));
    if max_amp > 0
        ylim([-max_amp*1.5, max_amp*1.5]);
    end
end

yyaxis right;
p2 = plot(t_axis, inst_freq_num, '-' ,'color', '#006400', 'DisplayName', '\textbf{Inst. Freq.}');
ylabel('\boldmath$\mathbf{Inst. \ Freq. \ [\omega_0]}$');

%green text
yline(0.6, '--', 'Color', '#006400', 'LineWidth', 3, 'HandleVisibility', 'off');
text(-23, 0.61, '\boldmath$\mathbf{\omega_{max}(SO)}$', 'Color', '#006400', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 14, 'Interpreter', 'latex');

ax = gca;
ax.YAxis(2).Color = 'k';

% Dynamic scaling for frequency
if ~isempty(ylim_plot1_freq)
    ylim(ylim_plot1_freq);
elseif any(mask_so)
    max_freq = max(abs(inst_freq_num(mask_so)));
    if max_freq > 0
        ylim([0, max_freq*1.5]);
    end
end

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
xlim([-27,27]);

legend('show', 'Location', 'best');

set(gca, 'XColor', 'k');
PlotUtils.styleAxes(gca);

% Align the 0 of both y-axes by minimally expanding limits
yyaxis left; yL = ylim;
yyaxis right; yR = ylim;
if yL(1) < 0 && yL(2) > 0 && yR(1) < 0 && yR(2) > 0
    rL = yL(2) / abs(yL(1));
    rR = yR(2) / abs(yR(1));
    Rs = [rL, rR, 1];
    best_R = Rs(1);
    min_penalty = inf;
    for R = Rs
        new_yL1 = -max(abs(yL(1)), yL(2)/R);
        new_yL2 = max(yL(2), abs(yL(1))*R);
        new_yR1 = -max(abs(yR(1)), yR(2)/R);
        new_yR2 = max(yR(2), abs(yR(1))*R);
        penalty = (new_yL2 - new_yL1)/(yL(2) - yL(1)) + (new_yR2 - new_yR1)/(yR(2) - yR(1));
        if penalty < min_penalty
            min_penalty = penalty;
            best_R = R;
        end
    end
    yyaxis left;  ylim([-max(abs(yL(1)), yL(2)/best_R), max(yL(2), abs(yL(1))*best_R)]);
    yyaxis right; ylim([-max(abs(yR(1)), yR(2)/best_R), max(yR(2), abs(yR(1))*best_R)]);
end

hold off;
