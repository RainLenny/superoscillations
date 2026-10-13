function fit_RLC()
    % FIT_RLC Fits theoretical RLC transfer function to experimental notch data.
    %
    % Circuit Model:
    %   Source Vin -> series R0 -> [N stages of series R_{k-1} and shunt Z_k] -> series RN -> RL (50 Ohm)
    %   where shunt branch k is a series RLC branch:
    %       Z_k(i*omega) = R_par_k + 1i*omega*L_k + 1/(1i*omega*C_k)
    %
    % Theoretical Transfer Function:
    %   H(i*omega) = Vout / Vin = (R_L / (R_N + R_L)) * prod_{k=1}^N (Z_k / (R_{k-1} + Z_k))
    %   Simulated spectrum analyzer power (dBm):
    %       P_sim(f) = 20 * log10(|H(2*pi*f)|) + P0

    %% 1. File Paths and Data Loading
    script_dir = fileparts(mfilename('fullpath'));
    data_dir   = fullfile(script_dir, 'data');
    mat_file   = fullfile(script_dir, 'optimized_RLC.mat');

    % --- Measured Component Values ---
    meas_file = fullfile(data_dir, 'measured_rlc_components.csv');
    opts_meas = detectImportOptions(meas_file);
    T_meas = readtable(meas_file, opts_meas);

    R_meas = T_meas{:, 2}.';         % Series resistors R_1 to R_N [Ohm]
    L_meas = T_meas{:, 3}.' * 1e-6;  % Inductances L_1 to L_N [H]
    C_meas = T_meas{:, 4}.' * 1e-12; % Capacitances C_1 to C_N [F]
    N = length(R_meas);
    R_L = 50; % Spectrum analyzer 50 Ohm load impedance

    % --- Measured Spectrum Analyzer Data ---
    sa_file = fullfile(data_dir, 'sa_notch.csv');
    try
        data_sa = readmatrix(sa_file, 'NumHeaderLines', 34);
    catch
        opts_sa = detectImportOptions(sa_file, 'NumHeaderLines', 34);
        data_tbl = readtable(sa_file, opts_sa);
        data_sa = [data_tbl{:, 1}, data_tbl{:, 2}];
    end

    valid_mask = ~isnan(data_sa(:, 1)) & ~isnan(data_sa(:, 2)) & (data_sa(:, 1) > 0);
    freq_full = data_sa(valid_mask, 1);
    amp_full  = data_sa(valid_mask, 2);

    %% 2. Frequency Range for Fitting
    fit_fmin = 200e3; % 200 kHz
    fit_fmax = 800e3; % 800 kHz

    fit_mask = (freq_full >= fit_fmin) & (freq_full <= fit_fmax);
    freq_fit = freq_full(fit_mask);
    amp_fit  = amp_full(fit_mask);

    %% 3. Parameter Bounds and Transformations
    % DMM tolerances around measured values:
    tol_R = 0.05; % +/- 5% for resistors
    tol_L = 0.10; % +/- 10% for inductors (including lead/core parasitics)
    tol_C = 0.10; % +/- 10% for capacitors (including stray capacitance)

    lb_R0   = 0.1;                      ub_R0   = 5.0;
    lb_R    = (1 - tol_R) * R_meas;     ub_R    = (1 + tol_R) * R_meas;
    lb_L    = (1 - tol_L) * L_meas;     ub_L    = (1 + tol_L) * L_meas;
    lb_C    = (1 - tol_C) * C_meas;     ub_C    = (1 + tol_C) * C_meas;
    lb_Rpar = 1.0 * ones(1, N);         ub_Rpar = 5000.0 * ones(1, N);

    lb = [lb_R0, lb_R, lb_L, lb_C, lb_Rpar];
    ub = [ub_R0, ub_R, ub_L, ub_C, ub_Rpar];

    % Smooth bounded mapping using sigmoid transformation: param = lb + (ub - lb)/(1 + exp(-u))
    decode = @(u) lb + (ub - lb) ./ (1 + exp(-u));
    encode = @(p) log((min(max(p, lb + 1e-4*(ub-lb)), ub - 1e-4*(ub-lb)) - lb) ./ ...
                      (ub - min(max(p, lb + 1e-4*(ub-lb)), ub - 1e-4*(ub-lb))));

    cost_fun = @(u) eval_cost(u, decode, freq_fit, amp_fit, R_L, N);

    %% 4. Multi-Start Optimization
    % Stage 1 corresponds to notch resonance near 540 kHz (R_par1 small ~ 2-5 Ohm)
    % Remaining stages are heavily damped (R_par ~ 500-2000 Ohm)
    rpar_guesses = {
        [2.5, 500  * ones(1, N-1)];
        [2.0, 1000 * ones(1, N-1)];
        [3.0, 2000 * ones(1, N-1)];
        [1.5, 300  * ones(1, N-1)]
    };

    fprintf('=======================================================\n');
    fprintf('  Fitting RLC Notch Filter (200 kHz - 800 kHz)\n');
    fprintf('=======================================================\n');
    fprintf('Number of stages N = %d, Load R_L = %.1f Ohm\n', N, R_L);

    best_cost = inf;
    best_u = [];

    opts = optimset('Display', 'none', 'MaxFunEvals', 25000, 'MaxIter', 10000);

    for g = 1:length(rpar_guesses)
        p_guess = [2.0, R_meas, L_meas, C_meas, rpar_guesses{g}];
        u0 = encode(p_guess);
        [u_opt, c_opt] = fminsearch(cost_fun, u0, opts);
        rmse_g = sqrt(c_opt);
        fprintf('  Candidate %d: RMSE = %.4f dB\n', g, rmse_g);
        if c_opt < best_cost
            best_cost = c_opt;
            best_u = u_opt;
        end
    end

    % Fine refinement pass
    opts_refine = optimset('Display', 'none', 'MaxFunEvals', 15000, 'MaxIter', 5000, 'TolFun', 1e-6, 'TolX', 1e-6);
    [best_u, best_cost] = fminsearch(cost_fun, best_u, opts_refine);

    %% 5. Unpack Optimized Parameters
    p_opt     = decode(best_u);
    R0_opt    = p_opt(1);
    R_opt     = p_opt(2:N+1);
    L_opt     = p_opt(N+2:2*N+1);
    C_opt     = p_opt(2*N+2:3*N+1);
    R_par_opt = p_opt(3*N+2:4*N+1);

    % Compute exact optimal baseline offset P0
    H_fit = calc_H(freq_fit, R0_opt, R_opt, L_opt, C_opt, R_par_opt, R_L);
    mag_fit_dB = 20 * log10(abs(H_fit));
    P0_final = mean(amp_fit(:) - mag_fit_dB(:));

    final_rmse = sqrt(best_cost);
    max_abs_err = max(abs((mag_fit_dB(:) + P0_final) - amp_fit(:)));

    fprintf('-------------------------------------------------------\n');
    fprintf('Optimization Succeeded:\n');
    fprintf('  Fit RMSE (200 - 800 kHz): %.4f dB\n', final_rmse);
    fprintf('  Max Absolute Error:       %.4f dB\n', max_abs_err);
    fprintf('  R0 = %.3f Ohm, P0 = %.3f dBm\n\n', R0_opt, P0_final);

    fprintf('Stage Parameters Summary:\n');
    fprintf(' Stage | R_opt [Ohm] | L_opt [uH] | C_opt [pF] | R_par [Ohm] | f_res [kHz]\n');
    fprintf('-------+-------------+------------+------------+-------------+------------\n');
    for k = 1:N
        f_res_k = 1 / (2 * pi * sqrt(L_opt(k) * C_opt(k))) / 1e3;
        fprintf('  %2d   |   %7.2f   |  %8.2f  |  %8.2f  |   %8.1f  |   %7.2f\n', ...
            k, R_opt(k), L_opt(k)*1e6, C_opt(k)*1e12, R_par_opt(k), f_res_k);
    end
    fprintf('-------------------------------------------------------\n');

    %% 6. Save Optimized Results
    H_opt = @(f) calc_H(f, R0_opt, R_opt, L_opt, C_opt, R_par_opt, R_L);
    save(mat_file, 'R0_opt', 'R_opt', 'L_opt', 'C_opt', 'R_par_opt', 'P0_final', 'H_opt', 'freq_full', 'R_L');
    fprintf('Saved optimized results to: %s\n', mat_file);

    %% 7. Plotting Fit Results
    H_full_sim = calc_H(freq_full, R0_opt, R_opt, L_opt, C_opt, R_par_opt, R_L);
    amp_full_sim = 20 * log10(abs(H_full_sim)) + P0_final;

    figure('Name', 'RLC Notch Filter Fit', 'Color', 'w');
    hold on;
    plot(freq_full / 1e3, amp_full, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Experiment (SA)');
    plot(freq_full / 1e3, amp_full_sim, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Theoretical Fit');

    xline(fit_fmin / 1e3, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Fit Bounds');
    xline(fit_fmax / 1e3, 'b:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    xlabel('Frequency [kHz]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude [dBm]', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('RLC Filter Fit [200-800 kHz] (RMSE = %.3f dB)', final_rmse), ...
        'FontSize', 13, 'FontWeight', 'bold');
    xlim([100, 950]);
    grid on;
    box on;
    legend('Location', 'southeast', 'FontSize', 11);
    hold off;
end

%% Helper: Least-Squares Cost Function
function cost = eval_cost(u, decode, freq_fit, amp_fit, R_L, N)
    p = decode(u);
    R0    = p(1);
    R     = p(2:N+1);
    L     = p(N+2:2*N+1);
    C     = p(2*N+2:3*N+1);
    R_par = p(3*N+2:4*N+1);

    H = calc_H(freq_fit, R0, R, L, C, R_par, R_L);
    mag_dB = 20 * log10(abs(H));
    P0 = mean(amp_fit(:) - mag_dB(:));
    residual = (mag_dB(:) + P0) - amp_fit(:);
    cost = mean(residual.^2);
end

%% Helper: Vectorized Transfer Function Calculation
function H = calc_H(f, R0, R, L, C, R_par, R_L)
    % CALC_H Vectorized calculation of the RLC ladder transfer function.
    %
    % Formula:
    %   H(i*omega) = (R_L / (R_N + R_L)) * prod_{k=1}^N (Z_k / (R_{k-1} + Z_k))
    %   where Z_k = R_par_k + 1i*omega*L_k + 1/(1i*omega*C_k)

    omega = 2 * pi * f(:).';              % 1 x M row vector of angular frequencies
    R_prev = [R0, R(1:end-1)].';          % N x 1 column vector of series resistors R_{k-1}

    % Shunt branch impedances Z_k (size: N x M via implicit expansion)
    Z = R_par(:) + 1i * omega .* L(:) + 1 ./ (1i * omega .* C(:));

    % Stage voltage dividers: Z_k / (R_{k-1} + Z_k)
    stage_H = Z ./ (R_prev + Z);

    % Total transfer function across all stages
    H_row = (R_L / (R(end) + R_L)) * prod(stage_H, 1);

    % Ensure output dimensions match the input frequency vector f
    H = reshape(H_row, size(f));
end
