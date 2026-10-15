%% fit_RLC.m
% Optimizes RLC notch filter parameters against experimental spectrum analyzer data.
% Enforces exact resonant frequency matching and removes resistor degeneracy.
clear; clc; close all;

%% 1. Paths & Data Loading
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(fileparts(script_dir))));
if exist('PlotUtils', 'class'), PlotUtils.setupDefaults(); end

data_file = fullfile(script_dir, 'data', 'sa_notch.csv');
mat_file  = fullfile(script_dir, 'optimized_RLC.mat');

raw_data  = readmatrix(data_file, 'NumHeaderLines', 34);
freq_full = raw_data(1:end-6, 1);
sa_dB     = raw_data(1:end-6, 2);

up_fit_bound = 700e3;
low_fit_bound = 200e3;
fit_mask  = (freq_full >= low_fit_bound) & (freq_full <= up_fit_bound);
f_fit     = freq_full(fit_mask);
dB_fit    = sa_dB(fit_mask);

%% 2. Exact Experimental Resonance & Capacitance Lock
[~, idx_min] = min(dB_fit);
win_idx      = max(1, idx_min-3) : min(numel(f_fit), idx_min+3);

if numel(win_idx) >= 3
    poly_c    = polyfit(f_fit(win_idx), dB_fit(win_idx), 2);
    f_res_exp = -poly_c(2) / (2 * poly_c(1));
    dB_min    = polyval(poly_c, f_res_exp);
else
    f_res_exp = f_fit(idx_min);
    dB_min    = dB_fit(idx_min);
end

calc_C = @(inductance) 1 ./ ((2 * pi * f_res_exp)^2 .* inductance);

%% 3. Fixed Parameters & Optimization Bounds for Active Variables [L, R_par, P0]
[R_1_fixed, R_2_fixed, R_L_fixed] = deal(50.0);
R_eq        = (R_1_fixed * (R_2_fixed + R_L_fixed)) / (R_1_fixed + R_2_fixed + R_L_fixed);
dB_pass_est = max(dB_fit);

R_par_init = R_eq / (10^((dB_pass_est - dB_min) / 20) - 1);
P0_init    = dB_pass_est - 20 * log10(R_L_fixed / (R_1_fixed + R_2_fixed + R_L_fixed));

vars_lb   = [80.0e-6,  2.0, 10.0];
vars_ub   = [120.0e-6, 30.0, 25.0];
vars_init = [100.0e-6, min(max(R_par_init, 5.0), 25.0), min(max(P0_init, 10.0), 25.0)];

%% 4. Circuit Model & Objective Function
% Full 7-parameter circuit vector: params_full = [R_1, R_2, L, C, R_par, P0, R_L]
calc_transfer_func = @(f, params_full) params_full(7) ./ (params_full(1) + params_full(2) + params_full(7) + ...
    params_full(1) * (params_full(2) + params_full(7)) ./ (params_full(5) + 1j*2*pi*f*params_full(3) + 1./(1j*2*pi*f*params_full(4))));

% Map the active optimization variables to the full 7-parameter physical circuit vector.
map_active_to_full = @(opt_vars) [R_1_fixed, R_2_fixed, opt_vars(1), calc_C(opt_vars(1)), opt_vars(2), opt_vars(3), R_L_fixed];

weights     = 1 + 4 * exp(-((f_fit - f_res_exp) / 40e3).^2);
residual_fn = @(opt_vars) 20 * log10(abs(calc_transfer_func(f_fit, map_active_to_full(opt_vars)))) + opt_vars(3) - dB_fit;

%% 5. Solve Optimization
if exist('lsqnonlin', 'file')
    map_scaled_to_active = @(s_scaled) vars_lb + s_scaled .* (vars_ub - vars_lb);
    opts = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
    s_scaled_opt = lsqnonlin(@(s_scaled) weights .* residual_fn(map_scaled_to_active(s_scaled)), (vars_init - vars_lb) ./ (vars_ub - vars_lb), zeros(1,3), ones(1,3), opts);
    vars_opt = map_scaled_to_active(s_scaled_opt);
else
    map_uncon_to_active = @(u_uncon) vars_lb + (vars_ub - vars_lb) ./ (1 + exp(-u_uncon));
    opts = optimset('Display', 'notify', 'MaxFunEvals', 10000, 'MaxIter', 5000, 'TolX', 1e-12, 'TolFun', 1e-12);
    u_uncon_opt = fminsearch(@(u_uncon) rms(weights .* residual_fn(map_uncon_to_active(u_uncon))), log((vars_init - vars_lb) ./ (vars_ub - vars_init)), opts);
    vars_opt = map_uncon_to_active(u_uncon_opt);
end

%% 6. Unpack & Save Results (Backward Compatible)
full_params_cell = num2cell(map_active_to_full(vars_opt));
[R_1_opt, R_2_opt, L_opt, C_opt, R_par_opt, P0_final, R_L_opt] = deal(full_params_cell{:});
[R0_opt, R_opt, R_L] = deal(R_1_opt, R_2_opt, R_L_opt);

f_res_opt = 1 / (2 * pi * sqrt(L_opt * C_opt));
rmse_opt  = rms(residual_fn(vars_opt));
H_opt     = @(f) calc_transfer_func(f, [full_params_cell{:}]);

save(mat_file, 'R_1_opt', 'R_2_opt', 'L_opt', 'C_opt', 'R_par_opt', 'P0_final', ...
               'H_opt', 'freq_full', 'R_L_opt', 'R0_opt', 'R_opt', 'R_L', 'f_res_opt');

%% 8. Verification Plot
figure('Name', 'RLC Notch Filter Fit');

hold on;
patch([low_fit_bound, up_fit_bound, up_fit_bound, low_fit_bound], [-10, -10, 15, 15], [0.9 0.95 1.0], ...
    'EdgeColor', 'none', 'DisplayName', 'Fit Range (200--800 kHz)');
plot(freq_full, sa_dB, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Experimental Data');
plot(freq_full, 20 * log10(abs(H_opt(freq_full))) + P0_final, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Fitted Model');
xline(f_res_opt, 'b:', 'LineWidth', 1.5, 'DisplayName', sprintf('Resonance (%.2f kHz)', f_res_opt / 1e3));

xlabel('Frequency [Hz]');
ylabel('Amplitude [dBm]');
xlim([0, 1e6]); ylim([-8, 10]);
legend('Location', 'southeast'); grid on;
if exist('PlotUtils', 'class'), PlotUtils.styleAxes(gca); end