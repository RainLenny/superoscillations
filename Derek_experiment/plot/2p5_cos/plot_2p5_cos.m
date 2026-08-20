clear; clc;
%% Importing signals
% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();
[VinSOFun_comp, angular_freqs, amps_SO] = generate_SO_from_dat_file('2p5_cos_Derek', 1, 1, false);


% Normalization such that the cos will be of value 1 at t=0
VinSOFun = @(t) real(VinSOFun_comp(t)) ./ real(VinSOFun_comp(0));


%% SAMPLING Constants
f_sampling = 600/(2*pi);
% fundamental_period of the superoscillating signal:
T_period = compute_fundamental_period([angular_freqs],f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample‐interval in seconds
t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1);
t_axis = t_axis(:);


%% Sample the signals
sampled_superoscillation = VinSOFun(t_axis);

%% Plot 1: SO in time with instantaneous frequency
% --- Plot 1 Y-Limit Controls ---
% Set to [min, max] for manual limits, or [] for dynamic scaling
ylim_plot1_so = [-4,7]*10^(-2)/0.0148;
ylim_plot1_freq = [-0.5,1.5];
% -------------------------------

figure;
[inst_freq_num, ~, ~] = compute_instantaneous_frequency(sampled_superoscillation, t_axis);

yyaxis left;
p1 = plot(t_axis, real(sampled_superoscillation), '-', 'Color', 'r', 'DisplayName', '\textbf{SO}');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');
ax = gca;
ax.YAxis(1).Color = 'k';

% Dynamic scaling for amplitude
mask_so = abs(t_axis) <= 5;
if ~isempty(ylim_plot1_so)
    ylim(ylim_plot1_so);
elseif any(mask_so)
    max_amp = max(abs(real(sampled_superoscillation(mask_so))));
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

%% Plot 2: SO zoom out

figure;

p1 = plot(t_axis, real(sampled_superoscillation), '-', 'Color', 'r', 'DisplayName', '\textbf{SO}');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');
ax = gca;
ax.YAxis(1).Color = 'k';




xlim([-70,70]);

legend('show', 'Location', 'best');

set(gca, 'XColor', 'k');
PlotUtils.styleAxes(gca);


xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');

hold off;



%% FFT the signals
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

% Calculate FFT and normalize by N for true amplitude
fft_superoscillation = fftshift(abs(fft(sampled_superoscillation,N))) / N;





%% Plot 2: SO FFT
figure;
hold on;
h1 = plot(freq_axis, fft_superoscillation, '-o', 'Color', 'r', 'DisplayName', '\textbf{SO}');

% 2. Vertical reference line & Annotation
xline(1.0, 'Color', 'black','LineWidth', 6);
text(1, 25, '\boldmath$\mathbf{\omega_0}$', ...
    'Color', 'k', 'FontSize', 18, 'Rotation', 90, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% Updated X and Y labels with bold math/text formatting
xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb.]}$');

legend([h1]);
grid off;
xlim([0,1.1])
% ylim([0,36])

set(gca, 'XColor', 'k');
PlotUtils.styleAxes(gca);

% % Reduce the number of y-ticks
% ax = gca;
% current_yticks = yticks;
% yticks(current_yticks(1:2:end));

hold off;