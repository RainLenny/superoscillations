%% PATH SETUP & CONSTANTS
clear; clc;
% Add all project subfolders to search path
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(project_root));

% Try to use PlotUtils if it exists, otherwise ignore
try
    PlotUtils.setupDefaults();
catch
end

%% 1. Define Signals
signal_scaling = 1.1;

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
try
    sig_configs = normalize_sig_configs(sig_configs, 'energy');
catch
    disp('Warning: normalize_sig_configs not found or failed, skipping normalization.');
end


%% 2. Run Simulations
% RLC parameters
omega0 = 1; % Resonant frequency
L = 1;
C0 = 1 / (omega0^2 * L); 
R = 0.02; % Small damping for exponential decay
alpha = 0.01; % Varactor nonlinearity (weakly nonlinear)

tspan = [0, 500];
y0 = [0; 0]; % initial charge and current

for i = 1:length(sig_configs)
    % The signal from generator is complex, we take the real part to drive the RLC circuit
    V_in = @(t) real(sig_configs(i).data(t));
    
    [tgrid, q, I, E_L, E_C, E_tot] = solve_RLC(L, R, C0, alpha, V_in, tspan, y0);
    
    sig_configs(i).tgrid = tgrid;
    sig_configs(i).E_tot = E_tot;
end

%% 3. Plot Results
figure;
ax1 = subplot(1,1,1);
hold on;

h_arr = gobjects(1, length(sig_configs));
for i = 1:length(sig_configs)
    h_arr(i) = plot(sig_configs(i).tgrid, sig_configs(i).E_tot, 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

grid on;
xlabel('Time [arb]');
ylabel('Total Energy $E_{total}(t)$', 'Interpreter', 'latex');
title('RLC Circuit Energy Dynamics');

try
    PlotUtils.styleAxes(ax1);
catch
end
legend(h_arr, 'Interpreter', 'latex');

%% Plot Signals
try
    plot_signals(sig_configs);
catch
end
