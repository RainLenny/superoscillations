%% Setup Parameters
clear; clc; close all;

%% 1. Global constants and Time
t_span = [0 600];

% Integration window parameters for J calculation
T1_integ = 235;   % fixed start time (T_start)

%% 2. Initial Conditions and System Constants
rho_init = [0; 0; -1];

params.T1         = 500;
params.T2         = 300;
params.rho30      = -1.0;
params.OmegaTilde = 0.01;

% Define the two TLS systems to compare
nu0_TLS1 = 1.00;
nu0_TLS2 = 1.01;

%% 3. Define Signals
signals_list = {};
signal_names = {};

% Generate raw signals (scaling arguments 1, 10 passed as before)
% Note: You no longer need to worry about manual scaling inside 'generate_signals' 
% because the normalization function will fix it.
[SO_signal, Cos_non_resosnant_signal, Cos_resosnant_signal] = generate_signals_Baranov_2014(1, 1);

signals_list{end+1} = SO_signal;
signal_names{end+1} = 'SO ';

signals_list{end+1} = Cos_non_resosnant_signal;
signal_names{end+1} = 'Cos 0.9[ω]';

% signals_list{end+1} = Cos_resosnant_signal;
% signal_names{end+1} = 'Cos 1[ω]';


% --- Apply Normalization ---
t_scan_for_norm = linspace(t_span(1), t_span(end), 5000);

signals_list = normalize_signal_list_by_energy(signals_list, t_scan_for_norm);
% --------------------------------

%% 4. Main Processing Loop
results = struct();

% Common time vector for interpolation
t_common = linspace(t_span(1), t_span(2), 2000)';

% Index corresponding to T1_integ
idx_start = find(t_common >= T1_integ, 1);
if isempty(idx_start), idx_start = 1; end

% Time axis for J plots (T2)
t_axis = t_common(idx_start:end);

for signal_idx = 1:length(signals_list)
    curr_signal = signals_list{signal_idx};
    disp(['Processing ' signal_names{signal_idx} '...']);

    % --- Run Solver for TLS 1 ---
    [t_out1, rho_out1] = solve_bloch_odes(t_span, rho_init, nu0_TLS1, params, curr_signal);

    % --- Run Solver for TLS 2 ---
    [t_out2, rho_out2] = solve_bloch_odes(t_span, rho_init, nu0_TLS2, params, curr_signal);

    % --- Interpolate to common grid ---
    rho1_interp = interp1(t_out1, rho_out1, t_common, 'pchip');
    rho2_interp = interp1(t_out2, rho_out2, t_common, 'pchip');

    % Store rho trajectories (for plotting later)
    results(signal_idx).rho_tls1 = rho1_interp;   % size: [length(t_common) x 3]
    results(signal_idx).rho_tls2 = rho2_interp;   % size: [length(t_common) x 3]
    results(signal_idx).name     = signal_names{signal_idx};

    % --- Calculate J for each component (rho1, rho2, rho3) ---
    J_vals = zeros(length(t_axis), 3); % columns: rho1, rho2, rho3

    for comp = 1:3
        r1 = rho1_interp(:, comp);
        r2 = rho2_interp(:, comp);

        numerator_integrand   = (r1 - r2).^2;
        denominator_integrand = 0.5 * (r1.^2 + r2.^2);

        num_int = cumtrapz(t_common(idx_start:end), numerator_integrand(idx_start:end));
        den_int = cumtrapz(t_common(idx_start:end), denominator_integrand(idx_start:end));

        den_int(den_int == 0) = eps;

        J_vals(:, comp) = num_int ./ den_int;
    end

    results(signal_idx).J = J_vals; % size: [length(t_axis) x 3]
end

%% Plot J_i separately (3 figures), LINEAR (not semilogy)
component_labels = {'\rho_1','\rho_2','\rho_3'};

for comp = 1:3
    figure('Color','w','Name',['J for ' component_labels{comp}]);
    hold on;

    for signal_idx = 1:length(signals_list)
        plot(t_axis, results(signal_idx).J(:, comp), 'LineWidth', 1.5, ...
            'DisplayName', results(signal_idx).name);
    end

    title(['Distinguishability J for ' component_labels{comp}]);
    xlabel('Integration Limit T_2');
    ylabel('J Parameter');
    grid on;
    legend('Location','best');
    xlim([T1_integ t_span(2)]);
end

%% Plot rho_i in separate figures (3 figures), each figure shows ALL signals
% Plots BOTH TLS1 and TLS2 for each signal (solid = TLS1, dashed = TLS2).
for comp = 1:3
    figure('Color','w','Name',['rho for ' component_labels{comp}]);
    hold on;

    for signal_idx = 1:length(signals_list)
        r_tls1 = results(signal_idx).rho_tls1(:, comp);
        r_tls2 = results(signal_idx).rho_tls2(:, comp);

        plot(t_common, r_tls1, 'LineWidth', 1.5, ...
            'DisplayName', [results(signal_idx).name ' (TLS1)']);
        plot(t_common, r_tls2, '--', 'LineWidth', 1.5, ...
            'DisplayName', [results(signal_idx).name ' (TLS2)']);
    end

    title(['Time-domain trajectories for ' component_labels{comp}]);
    xlabel('Time t');
    ylabel(['Component ' component_labels{comp}]);
    grid on;
    legend('Location','best');
    xlim(t_span);
end

%% Plot all input signals
figure('Color','w','Name','Input Signals');
hold on;

for signal_idx = 1:length(signals_list)
    curr_signal = signals_list{signal_idx};
    f_vals = arrayfun(curr_signal, t_common);
    plot(t_common, f_vals, 'LineWidth', 1.5, ...
        'DisplayName', signal_names{signal_idx});
end

xlabel('Time t');
ylabel('Signal amplitude');
title('Comparison of all input signals');
grid on;
legend('Location','best');
xlim(t_span);
