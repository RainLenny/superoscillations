%% fit_RLC.m
% Fit RLC notch filter parameters to experimental spectrum analyzer data (sa_notch.csv).
% Optimization window: 200 kHz - 800 kHz.
%
% Circuit topology:
%   Z_branch = i*omega*L + 1/(i*omega*C) + R_par
%   Z = Z_branch || R_L = (Z_branch * R_L) / (Z_branch + R_L)
%   H(i*omega) = V_out / V_in = Z / (R + Z)
%
% Output:
%   Generates 'optimized_RLC.mat' containing:
%   'R0_opt', 'R_opt', 'L_opt', 'C_opt', 'R_par_opt', 'P0_final', 'H_opt', 'freq_full', 'R_L'

clear; clc;

%% 1. Path Setup & Load Experimental Data
script_dir = fileparts(mfilename('fullpath'));
data_file  = fullfile(script_dir, 'data', 'sa_notch.csv');
mat_file   = fullfile(script_dir, 'optimized_RLC.mat');

% Read spectrum analyzer data (skip 34 header lines, remove trailing metadata)
raw_data = readmatrix(data_file, 'NumHeaderLines', 34);
valid_idx = ~isnan(raw_data(:, 1)) & ~isnan(raw_data(:, 2));
freq_full = raw_data(valid_idx, 1);
dB_full   = raw_data(valid_idx, 2);

%% 2. Define Fitting Range (200 kHz to 800 kHz)
fit_mask = (freq_full >= 200e3) & (freq_full <= 800e3);
f_fit    = freq_full(fit_mask);
dB_fit   = dB_full(fit_mask);

%% 3. Nominal Circuit Components & Physical Bounds
% Given component specifications:
R0        = 46.4;     % Resistor [Ohm]
L0        = 100e-6;   % Inductor [H] (100 uH)
C0        = 870e-12;  % Capacitor [F] (870 pF)
R_L0      = 50;       % Load resistor [Ohm]
R_par_max = 20;       % Max parasitic resistance [Ohm]

% Component bounds:
lb = [0.8 * R0, 0.8 * L0, 0.8 * C0, 0, 0.8 * R_L0];
ub = [500 * R0, 1.2 * L0, 1.2 * C0, R_par_max, 1.2 * R_L0];

%% 4. Circuit Transfer Function Formulation
% Z_branch = i*omega*L + 1/(i*omega*C) + R_par
% Z = Z_branch || R_L
% H(i*omega) = Z / (R + Z)
calc_H = @(f, R, L, C, R_par, RL) ...
    (( (1j*2*pi*f*L + 1./(1j*2*pi*f*C) + R_par) .* RL ) ./ ...
     ( (1j*2*pi*f*L + 1./(1j*2*pi*f*C) + R_par) + RL )) ./ ...
    ( R + (( (1j*2*pi*f*L + 1./(1j*2*pi*f*C) + R_par) .* RL ) ./ ...
          ( (1j*2*pi*f*L + 1./(1j*2*pi*f*C) + R_par) + RL )) );

%% 5. Bounded Optimization (Nelder-Mead with Variable Transformation)
% Smooth trigonometric parameter mapping: p = lb + (ub - lb) * (sin(theta) + 1)/2
to_param = @(theta) lb + (ub - lb) .* (sin(theta) + 1) / 2;
to_theta = @(p) asin(min(max(2 * (p - lb) ./ (ub - lb) - 1, -1), 1));

% Initial guess (centered in allowable range)
p_init     = [R0, L0, C0, 2.5, R_L0];
theta_init = to_theta(p_init);

% Objective function: RMSE in dB with analytical baseline offset P0
cost_fun = @(theta) norm( ...
    ( 20*log10(abs(calc_H(f_fit, ...
        lb(1) + (ub(1)-lb(1))*(sin(theta(1))+1)/2, ...
        lb(2) + (ub(2)-lb(2))*(sin(theta(2))+1)/2, ...
        lb(3) + (ub(3)-lb(3))*(sin(theta(3))+1)/2, ...
        lb(4) + (ub(4)-lb(4))*(sin(theta(4))+1)/2, ...
        lb(5) + (ub(5)-lb(5))*(sin(theta(5))+1)/2))) + ...
      mean(dB_fit - 20*log10(abs(calc_H(f_fit, ...
        lb(1) + (ub(1)-lb(1))*(sin(theta(1))+1)/2, ...
        lb(2) + (ub(2)-lb(2))*(sin(theta(2))+1)/2, ...
        lb(3) + (ub(3)-lb(3))*(sin(theta(3))+1)/2, ...
        lb(4) + (ub(4)-lb(4))*(sin(theta(4))+1)/2, ...
        lb(5) + (ub(5)-lb(5))*(sin(theta(5))+1)/2)))) ) - dB_fit );

opts = optimset('Display', 'notify', 'TolX', 1e-6, 'TolFun', 1e-6, ...
                'MaxFunEvals', 40000, 'MaxIter', 40000);
theta_opt = fminsearch(cost_fun, theta_init, opts);
p_opt     = to_param(theta_opt);

%% 6. Extract Optimized Parameters
R_opt     = p_opt(1);
R0_opt    = R_opt;           % Canonical alias for series resistor
L_opt     = p_opt(2);
C_opt     = p_opt(3);
R_par_opt = p_opt(4);
R_L       = p_opt(5);        % Optimized load resistor
f_res_opt = 1 / (2 * pi * sqrt(L_opt * C_opt));

% Analytical optimal baseline offset P0 [dBm]
P0_final  = mean(dB_fit - 20*log10(abs(calc_H(f_fit, R_opt, L_opt, C_opt, R_par_opt, R_L))));

%% 7. Construct Self-Contained Transfer Function Handle
% Constructed using explicit numerical values so that it loads cleanly without external m-file dependency
H_opt = @(f) (((1j*2*pi*f*L_opt + 1./(1j*2*pi*f*C_opt) + R_par_opt) * R_L) ./ ...
    ((1j*2*pi*f*L_opt + 1./(1j*2*pi*f*C_opt) + R_par_opt) + R_L)) ./ ...
    (R_opt + (((1j*2*pi*f*L_opt + 1./(1j*2*pi*f*C_opt) + R_par_opt) * R_L) ./ ...
    ((1j*2*pi*f*L_opt + 1./(1j*2*pi*f*C_opt) + R_par_opt) + R_L)));

%% 8. Save Optimized Results to MAT-file
save(mat_file, 'R0_opt', 'R_opt', 'L_opt', 'C_opt', 'R_par_opt', 'P0_final', 'H_opt', 'freq_full', 'R_L', 'f_res_opt');
fprintf('Successfully saved optimized parameters to:\n  %s\n\n', mat_file);

%% 9. Display Fit Summary
fprintf('================ Optimization Summary ================\n');
fprintf('  R_opt      = %8.4f Ohm      (nominal: %5.1f Ohm)\n', R_opt, R0);
fprintf('  L_opt      = %8.4e H (%.2f uH, nominal: %5.1f uH)\n', L_opt, L_opt*1e6, L0*1e6);
fprintf('  C_opt      = %8.4e F (%.2f pF, nominal: %5.1f pF)\n', C_opt, C_opt*1e12, C0*1e12);
fprintf('  R_par_opt  = %8.4f Ohm      (max bound: %5.1f Ohm)\n', R_par_opt, R_par_max);
fprintf('  R_L        = %8.4f Ohm      (nominal: %5.1f Ohm)\n', R_L, R_L0);
fprintf('  f_res_opt  = %8.2f kHz\n', f_res_opt / 1e3);
fprintf('  P0_final   = %8.4f dBm\n', P0_final);
fit_rmse = norm((20*log10(abs(H_opt(f_fit))) + P0_final) - dB_fit) / sqrt(length(f_fit));
fprintf('  Fit RMSE   = %8.4f dB   (over 200 kHz - 800 kHz)\n', fit_rmse);
fprintf('======================================================\n');
