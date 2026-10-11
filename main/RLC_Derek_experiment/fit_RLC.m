function fit_RLC()
    % Read data
    filename = 'Derek_experiment/data/sa_notch.csv';
    
    % Use readmatrix to automatically handle header strings and get numerical data
    % Data starts from line 35, we can just read the whole file and grab the bottom part
    % Or use readtable with 'NumHeaderLines'
    try
        opts = detectImportOptions(filename, 'NumHeaderLines', 34);
        data = readtable(filename, opts);
        freq = data.Var1;
        amp_meas = data.Var2;
    catch
        % Fallback for older MATLAB versions
        data = readmatrix(filename, 'NumHeaderLines', 34);
        freq = data(:,1);
        amp_meas = data(:,2);
    end
    
    % Full data arrays for plotting (ignore DC spike and NaNs)
    valid_full = (freq > 0) & ~isnan(freq) & ~isnan(amp_meas);
    freq_full = freq(valid_full);
    amp_full = amp_meas(valid_full);
    
    % --- FIT RANGE SETTINGS ---
    % Define the frequency range over which to optimize (in Hz)
    fit_fmin = 2e5;       % Minimum frequency for fit
    fit_fmax = 9e5;    % Maximum frequency for fit
    
    % Extract the subset used for the optimizer
    valid_fit = (freq_full >= fit_fmin) & (freq_full <= fit_fmax);
    freq_fit = freq_full(valid_fit);
    amp_fit = amp_full(valid_fit);
    
    N = 9;
    R_L = 50; % Load resistor (spectrum analyzer impedance)
    
    % Find notch frequency roughly (within fit range)
    [~, min_idx] = min(amp_fit);
    f0_guess = freq_fit(min_idx);
    
    % Step 1: Global fit for nominal values
    % Params: [Z0, R_par_nom, P0, f0]
    
    function amp_sim = sim_uniform(p, f)
        Z0 = p(1);
        R_par = p(2);
        P0 = p(3);
        f0 = p(4);
        
        L = Z0 / (2*pi*f0);
        C = 1 / (2*pi*f0*Z0);
        
        R_vec = 47 * ones(1, N);
        L_vec = L * ones(1, N);
        C_vec = C * ones(1, N);
        R_par_vec = R_par * ones(1, N);
        
        H = calc_H(f, R_vec, L_vec, C_vec, R_par_vec, R_L);
        amp_sim = 20*log10(abs(H)) + P0;
    end

    function cost = cost_uniform(p)
        amp_sim = sim_uniform(p, freq_fit);
        cost = sum((amp_sim(:) - amp_fit(:)).^2);
    end

    % Global optimization for Z0
    Z0_test = logspace(0, 3, 50);
    best_cost = inf;
    best_Z0 = 50;
    for z = Z0_test
        p_test = [z, 1, 26, f0_guess]; 
        c = cost_uniform(p_test);
        if c < best_cost
            best_cost = c;
            best_Z0 = z;
        end
    end
    
    p0_global = [best_Z0, 1.0, 26, f0_guess];
    lb_global = [best_Z0/10, 0.01, 0, f0_guess*0.5];
    ub_global = [best_Z0*10, 100, 50, f0_guess*1.5];
    
    options = optimoptions('fmincon','Display','none');
    p_opt_global = fmincon(@cost_uniform, p0_global, [], [], [], [], lb_global, ub_global, [], options);
    
    Z0_opt = p_opt_global(1);
    R_par_nom = p_opt_global(2);
    P0_opt = p_opt_global(3);
    f0_opt = p_opt_global(4);
    
    L_nom = Z0_opt / (2*pi*f0_opt);
    C_nom = 1 / (2*pi*f0_opt*Z0_opt);
    
    fprintf('Nominal values found:\n');
    fprintf('L = %g H\n', L_nom);
    fprintf('C = %g F\n', C_nom);
    fprintf('R_par = %g Ohm\n', R_par_nom);
    fprintf('P0 = %g dB\n', P0_opt);
    fprintf('f0 = %g Hz\n', f0_opt);
    
    % Step 2: Full optimization
    % To make optimization robust, we optimize normalized parameters x near 1
    % x = [r1..r9, l1..l9, c1..c9, rpar1..rpar9, p0_scaled]
    % r_k = R_k / 47
    % l_k = L_k / L_nom
    % c_k = C_k / C_nom
    % rpar_k = R_par_k / R_par_nom
    % p0_scaled = P0 / P0_opt
    
    x0 = ones(1, 4*N + 1);
    
    % Allow 20% tolerance for R, 20% for L and C, and bounds for Rpar and P0
    lb = [0.8*ones(1,N), 0.8*ones(1,N), 0.8*ones(1,N), 0.01*ones(1,N), 0.5];
    ub = [1.2*ones(1,N), 1.2*ones(1,N), 1.2*ones(1,N), 100*ones(1,N), 1.5];
    
    function amp_sim = sim_full(x, f)
        R_vec = x(1:N) * 47;
        L_vec = x(N+1:2*N) * L_nom;
        C_vec = x(2*N+1:3*N) * C_nom;
        R_par_vec = x(3*N+1:4*N) * R_par_nom;
        P0 = x(end) * P0_opt;
        
        H = calc_H(f, R_vec, L_vec, C_vec, R_par_vec, R_L);
        amp_sim = 20*log10(abs(H)) + P0;
    end

    function cost = cost_full(x)
        amp_sim = sim_full(x, freq_fit);
        cost = sum((amp_sim(:) - amp_fit(:)).^2);
    end
    
    fprintf('\nStarting full optimization...\n');
    options_full = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 5000);
    x_opt = fmincon(@cost_full, x0, [], [], [], [], lb, ub, [], options_full);
    
    R_opt = x_opt(1:N) * 47;
    L_opt = x_opt(N+1:2*N) * L_nom;
    C_opt = x_opt(2*N+1:3*N) * C_nom;
    R_par_opt = x_opt(3*N+1:4*N) * R_par_nom;
    P0_final = x_opt(end) * P0_opt;
    
    fprintf('\nOptimized Values:\n');
    disp('R (Ohms):'); disp(R_opt);
    disp('L (Henries):'); disp(L_opt);
    disp('C (Farads):'); disp(C_opt);
    disp('R_par (Ohms):'); disp(R_par_opt);
    fprintf('P0 = %g dB\n', P0_final);
    
    % Define optimized transfer function as a function of w
    H_opt = @(w) calc_H(w, R_opt, L_opt, C_opt, R_par_opt, R_L);
    
    % Save values to file
    save('optimized_RLC.mat', 'R_opt', 'L_opt', 'C_opt', 'R_par_opt', 'P0_final', 'H_opt', 'freq_full');
    
    % Plot results
    figure;
    plot(freq_full, amp_full, 'k-', 'LineWidth', 1.5); hold on;
    plot(freq_full, sim_uniform(p_opt_global, freq_full), 'b--', 'LineWidth', 1.5);
    plot(freq_full, sim_full(x_opt, freq_full), 'r-', 'LineWidth', 1.5);
    % Optional: Highlight the region used for fitting
    xline(fit_fmin, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
    xline(fit_fmax, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (dBm)');
    legend('Measurement', 'Nominal Fit', 'Full Fit');
    title('RLC Notch Filter Fit');
    grid on;
end

function H = calc_H(f, R, L, C, R_par, R_L)
    N = length(R);
    omega = 2 * pi * f(:).'; % Force row vector to avoid matrix expansion

    
    % Vectorize over frequencies
    M = length(omega);
    Z = zeros(N, M);
    
    for k = N:-1:1
        % Admittance of the shunt branch
        Y_shunt = 1 ./ (1i * omega * L(k) + 1 ./ (1i * omega * C(k)) + R_par(k));
        
        if k == N
            % Z_N = Z_out
            Z(k, :) = 1 ./ (1 / R_L + Y_shunt);
        else
            Z(k, :) = 1 ./ (Y_shunt + 1 ./ (R(k+1) + Z(k+1, :)));
        end
    end
    
    H = Z(1, :) ./ (R(1) + Z(1, :));
    for k = 2:N
        H = H .* (Z(k, :) ./ (R(k) + Z(k, :)));
    end
end
