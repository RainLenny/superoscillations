%% PATH SETUP & CONSTANTS
clear; clc;
% Add all project subfolders to search path
addpath(genpath(fileparts(mfilename('fullpath'))));

PlotUtils.setupDefaults();

[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov');
[Cos_signal_unscaled, angular_freqs_COS] = generate_Cos_reference(0.7, -1, true);

% Scale Cos_signal symbolically to match the peak of SO_signal
[~, peak_SO] = normalize_signals(SO_signal, 'peak');
[~, peak_COS] = normalize_signals(Cos_signal_unscaled, 'peak');
scale_factor = peak_SO / peak_COS;

[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.7, -scale_factor, true);

% Normalize both signals symbolically by the peak of the first signal
[SO_signal, Cos_signal] = normalize_signals({SO_signal, Cos_signal}, 'peak');

nu0 = 1;
tspan = [0,500];
Jtot = 0.03;


[tgrid_SO, Pe_SO] = JC_drive_only(nu0, Jtot, SO_signal, tspan);
[tgrid_COS, Pe_COS] = JC_drive_only(nu0, Jtot, Cos_signal, tspan);

%% PLOTS
figure;
ax1 = nexttile;
hold on;

plot(tgrid_SO, Pe_SO, 'LineWidth',5,'Color','red');
plot(tgrid_COS, Pe_COS, 'LineWidth',5,'Color','blue');

grid on;
xlabel('Time [arb]');
ylabel('Excited state probability');
PlotUtils.styleAxes(ax1);
legend({'SO','COS'});


plot_signals(SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS)