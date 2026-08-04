%% PATH SETUP & CONSTANTS
clear; clc;
% Add all project subfolders to search path
addpath(genpath(fileparts(mfilename('fullpath'))));

PlotUtils.setupDefaults();

%% 1. Define Signals
signal_scaling = 1;

sig_configs = struct('name', {}, 'data', {}, 'color', {}, 'freqs', {});

% --- Signal 1: SO ---
[SO_signal, angular_freqs_SO] = generate_SO_from_dat_file('SO_Baranov', 1, signal_scaling, true, true);
sig_configs(1).name = '\textbf{SO }';
sig_configs(1).data = SO_signal;
sig_configs(1).color = 'r';
sig_configs(1).freqs = angular_freqs_SO;

% --- Signal 2: COS ---
[Cos_signal, angular_freqs_COS] = generate_Cos_reference(0.9, signal_scaling, true, true);
sig_configs(2).name = '\boldmath$\mathbf{0.9\omega_0}$';
sig_configs(2).data = Cos_signal;
sig_configs(2).color = 'b';
sig_configs(2).freqs = angular_freqs_COS;

% Normalize all signals symbolically by the peak of the first signal
sig_configs = normalize_sig_configs(sig_configs, 'energy');


%% 2. Run Simulations
nu0 = 1;
tspan = [0,500];
Jtot = 0.02;

for i = 1:length(sig_configs)
    [tgrid, Pe] = JC_drive_only(nu0, Jtot, sig_configs(i).data, tspan);
    sig_configs(i).tgrid = tgrid;
    sig_configs(i).Pe = Pe;
end

%% 3. Plot Results
figure;
ax1 = nexttile;
hold on;

h_arr = gobjects(1, length(sig_configs));
for i = 1:length(sig_configs)
    h_arr(i) = plot(sig_configs(i).tgrid, sig_configs(i).Pe, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

grid on;
xlabel('Time [arb]');
ylabel('Excited state probability');
PlotUtils.styleAxes(ax1);
legend(h_arr);

%% Plot Signals
plot_signals(sig_configs);