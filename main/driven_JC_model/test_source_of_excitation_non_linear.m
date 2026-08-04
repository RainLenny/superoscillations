%% PATH SETUP & CONSTANTS
clear; clc; close all;


%% SETTINGS FOR PARAMETER SWEEP
Jtots = logspace(log10(0.005), log10(0.025), 100);

%% CONSTANTS

% Time and envelope parameters
t_0 = 250; 
T   = 100;

% Physical and simulation parameters
nu0      = 1.00;
Ce0      = 0; % Initial conditions of the TLS (ground state)


T_final = 500;

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


% Normalize all signals symbolically by the energy of the first signal
sig_configs = normalize_sig_configs(sig_configs, 'energy');

%% Dynamics computation
tspan = [0, T_final];


for i = 1:length(sig_configs)
    sig_configs(i).sigma_z_max = zeros(size(Jtots));
    sig_configs(i).sigma_x_max = zeros(size(Jtots));
    sig_configs(i).sigma_y_max = zeros(size(Jtots));
end

for j = 1:length(Jtots)
    Jtot = Jtots(j);
    
    for i = 1:length(sig_configs)
        
        [~, Pe, Ce, Cg] = JC_drive_only(nu0, Jtot, sig_configs(i).data, tspan, Ce0);
        
        % Calculate observables
        % sigma_z = |Ce|^2 - |Cg|^2 = 2Pe - 1
        % sigma_x = 2 * Re(Cg^* Ce)
        % sigma_y = 2 * Im(Cg^* Ce)
        
        sigma_z = 2 * Pe - 1;
        sigma_x = 2 * real(conj(Cg) .* Ce);
        sigma_y = 2 * imag(conj(Cg) .* Ce);
        
        % Record the maximum of the observables
        sig_configs(i).sigma_z_max(j) = max(Pe); % Proxy for max change in sigma_z
        sig_configs(i).sigma_x_max(j) = max(abs(sigma_x));
        sig_configs(i).sigma_y_max(j) = max(abs(sigma_y));
    end

    fprintf('Completed Jtot = %.2f\n', Jtot);
end


%% PLOTS
PlotUtils.setupDefaults();

% Set default interpreters to LaTeX to ensure custom formatting functions flawlessly
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% --- Figure 1: Amplitude Scaling Log-Log Plot (sigma_z proxy, max Pe) ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test (max P_e)');
hold on;

for i = 1:length(sig_configs)
    loglog(Jtots, sig_configs(i).sigma_z_max, '-', 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

% Plot reference slopes
ref_Jtot = Jtots;
% Align references to the first index of SO (first signal) for visualization
ref_1st_order = sig_configs(1).sigma_z_max(1) * (ref_Jtot / ref_Jtot(1)).^2;

loglog(ref_Jtot, ref_1st_order, '--k', 'DisplayName', '\boldmath$\propto J_{tot}^2$ \textbf{(1st order)}', 'LineWidth', 5);

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel('\boldmath$\max(P_e)$ \textbf{(Max Excitation Probability)}');

ylim([-inf, 1])

legend('show', 'Location', 'best');

% Apply Bold LaTeX Axis Ticks
PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold'); % Enforce bold ticks/scales post-styling


% --- Figure 2: Log-Log Slopes (max Pe) ---
figure('Color', 'w', 'Name', 'Log-Log Slopes (max P_e)');
hold on;

Jtot_mid = exp((log(Jtots(1:end-1)) + log(Jtots(2:end))) / 2);

for i = 1:length(sig_configs)
    slope_val = diff(log(sig_configs(i).sigma_z_max)) ./ diff(log(Jtots));
    semilogx(Jtot_mid, slope_val, '-', 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

yline(2, '--k', 'DisplayName', '\textbf{1st order (slope=2)}', 'LineWidth', 5);

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel('\boldmath$d(\log P_e) / d(\log J_{tot})$ \textbf{(Slope)}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'FontWeight', 'bold');


% --- Figure 3: Amplitude Scaling Log-Log Plot (sigma_x) ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test (sigma_x)');
hold on;

for i = 1:length(sig_configs)
    loglog(Jtots, sig_configs(i).sigma_x_max, '-', 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

% Plot reference slopes
ref_1st_order_rho = sig_configs(1).sigma_x_max(1) * (ref_Jtot / ref_Jtot(1)).^1;

loglog(ref_Jtot, ref_1st_order_rho, '--k', 'DisplayName', '\boldmath$\propto J_{tot}^1$ \textbf{(1st order)}', 'LineWidth', 5);

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel('\boldmath$\max(|\langle\sigma_x\rangle|)$');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold');


% --- Figure 4: Log-Log Slopes (sigma_x) ---
figure('Color', 'w', 'Name', 'Log-Log Slopes (sigma_x)');
hold on;

for i = 1:length(sig_configs)
    slope_val = diff(log(sig_configs(i).sigma_x_max)) ./ diff(log(Jtots));
    semilogx(Jtot_mid, slope_val, '-', 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

yline(1, '--k', 'DisplayName', '\textbf{1st order (slope=1)}', 'LineWidth', 5);

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel('\boldmath$d(\log |\langle\sigma_x\rangle|) / d(\log J_{tot})$ \textbf{(Slope)}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'FontWeight', 'bold');


% --- Figure 5: Amplitude Scaling Log-Log Plot (sigma_y) ---
figure('Color', 'w', 'Name', 'Amplitude Scaling Test (sigma_y)');
hold on;

for i = 1:length(sig_configs)
    loglog(Jtots, sig_configs(i).sigma_y_max, '-', 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

% Plot reference slopes
ref_1st_order_rho2 = sig_configs(1).sigma_y_max(1) * (ref_Jtot / ref_Jtot(1)).^1;

loglog(ref_Jtot, ref_1st_order_rho2, '--k', 'DisplayName', '\boldmath$\propto J_{tot}^1$ \textbf{(1st order)}', 'LineWidth', 5);

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel('\boldmath$\max(|\langle\sigma_y\rangle|)$');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold');


% --- Figure 6: Log-Log Slopes (sigma_y) ---
figure('Color', 'w', 'Name', 'Log-Log Slopes (sigma_y)');
hold on;

for i = 1:length(sig_configs)
    slope_val = diff(log(sig_configs(i).sigma_y_max)) ./ diff(log(Jtots));
    semilogx(Jtot_mid, slope_val, '-', 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end

yline(1, '--k', 'DisplayName', '\textbf{1st order (slope=1)}', 'LineWidth', 5);

grid on;
xlabel(' \textbf{Coupling (signal magnitude)}');
ylabel('\boldmath$d(\log |\langle\sigma_y\rangle|) / d(\log J_{tot})$ \textbf{(Slope)}');
legend('show', 'Location', 'best');

PlotUtils.styleAxes(gca);
set(gca, 'XScale', 'log', 'FontWeight', 'bold');
