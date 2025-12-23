function [tgrid, Pe] = old_multimode_JC_driven(omega, nu0, nmax, J0, Drive_integral, T_final)
%MULTIMODE_JC_DRIVEN Time-evolve a driven multi-mode Jaynes–Cummings model
%
%   [tgrid, Pe] = multimode_JC_driven(omega, nu0, nmax, J0, Drive_integral, T_final)
%
%   INPUTS:
%     omega           : 1×K vector of mode frequencies (or similar), used via Delta = nu0 - omega
%     nu0             : cavity / mode frequency (scalar)
%     nmax            : maximum photon number per mode
%     J0              : coupling strength (scalar)
%     Drive_integral  : function handle of the drive integral, Drive_integral(t)
%     T_final         : final simulation time
%
%   OUTPUTS:
%     tgrid : time points from ODE solver
%     Pe    : excited-state population of the TLS at each time in tgrid

    %% -------------------- NUMERICAL PARAMETERS -------------------------
    RelTol = 1e-9;            % relative tolerance for ODE solver
    AbsTol = 1e-9;            % absolute tolerance for ODE solver
    tspan  = [0, T_final];    % time interval

    %% --------------------- PHYSICAL PARAMETERS -------------------------
    Delta = nu0 - omega;      % detuning for each mode
    K     = numel(Delta);     % number of modes

    %% ------------------------- DIMENSIONS ------------------------------
    dim_ph1 = nmax + 1;       % local dimension per mode (0..nmax photons)
    dim_ph  = dim_ph1^K;      % total photonic Hilbert space dimension
    dim_tls = 2;              % two-level system dimension
    dim_tot = dim_ph * dim_tls;

    fprintf('K=%d modes, nmax=%d -> dim_ph=%d, dim_tot=%d\n', ...
            K, nmax, dim_ph, dim_tot);

    %% ------------------------- TLS OPERATORS ---------------------------
    % Raising/lowering operators in TLS space
    sp     = sparse([0 0; 1 0]);   % |e><g|
    sm     = sp.';                 % |g><e|
    Id_tls = speye(dim_tls);
    Id_ph  = speye(dim_ph);

    % Extend TLS operators to full Hilbert space (photons ⊗ TLS)
    Splus  = kron(Id_ph, sp);
    Sminus = kron(Id_ph, sm);

    %% -------------------- MODE ANNIHILATION OPERATORS -----------------
    % A_ops{k} acts as annihilation operator on mode k ⊗ identity on others ⊗ TLS
    A_ops = multimode_annihilation_ops_sparse(nmax, K, dim_tls);

    % Number operators for each mode: N_k = a_k^† a_k
    N_ops = cell(K,1);
    for k = 1:K
        Ak = A_ops{k};
        N_ops{k} = Ak' * Ak;
    end

    %% --------------------- STATIC HAMILTONIAN H0 ----------------------
    % H0 = sum_k [ -Delta_k * N_k + J0 * (a_k S^+ + a_k^† S^-) ]
    H0 = sparse(dim_tot, dim_tot);
    for k = 1:K
        Ak  = A_ops{k};
        Akd = Ak';
        H0  = H0 ...
              + (-Delta(k)) * N_ops{k} ...
              + J0 * (Ak * Splus + Akd * Sminus);
    end

    %% ------------------------ INITIAL STATE ----------------------------
    % Photonic vacuum in all modes
    vac_ph = sparse(dim_ph, 1);
    vac_ph(1) = 1;                % |0,0,...,0>

    % TLS ground state |g> = [1; 0]
    g_tls  = sparse([1; 0]);

    % Total initial state ψ(0) = |vac_ph> ⊗ |g>
    psi0   = kron(vac_ph, g_tls);

    %% ---------------- INDEX SET FOR EXCITED-STATE POPULATION ----------
    % TLS basis ordering assumed: [|g>; |e>] per photonic state,
    % so all excited components are entries 2, 4, 6, ...
    idx_e = 2:2:dim_tot;

    %% ---------------------- DRIVE FUNCTION f(t) ------------------------
    % Time-dependent drive coefficient coupling to S^+ and S^-
    f = @(t) -1i * J0 * Drive_integral(t) .* exp(1i * nu0 * t);

    %% --------------------- SCHRÖDINGER RHS HANDLE ---------------------
    % dψ/dt = -i [ H0 + f(t) S^+ + conj(f(t)) S^- ] ψ
    ode_rhs = @(t, psi) schrodinger_rhs_precomputed(t, psi, H0, Splus, Sminus, f);

    %% ---------------------- SPARSITY PATTERN --------------------------
    % Pattern used by ODE solver for internal Jacobian
    Jpattern = spones(H0 + Splus + Sminus);

    %% --------------------- OUTPUT FUNCTION (PROGRESS) -----------------
    output_fun = @(t, y, flag) local_output_fun(t, y, flag, T_final);

    %% -------------------------- ODE OPTIONS ---------------------------
    opts = odeset('RelTol',  RelTol, ...
                  'AbsTol',  AbsTol, ...
                  'JPattern', Jpattern, ...
                  'OutputFcn', output_fun);

    %% --------------------------- ODE SOLVE ----------------------------
    [tgrid, psi_all] = ode15s(ode_rhs, tspan, psi0, opts);

    %% -------------------- EXCITED-STATE POPULATION --------------------
    % Pe(t) = sum over all photonic configs of |⟨e, n | ψ(t)⟩|^2
    Pe = sum(abs(psi_all(:, idx_e)).^2, 2);

    fprintf('Progress: 100%%\n');
end

%% ======================== LOCAL FUNCTIONS ============================

function status = local_output_fun(t_curr, ~, flag, T_final)
%LOCAL_OUTPUT_FUN Progress printing for ODE solver

    persistent last_frac

    switch flag
        case 'init'
            last_frac = -inf;
            status = 0;
            return

        case ''  % during integration steps
            if isempty(t_curr)
                status = 0;
                return
            end

            t_now = t_curr(end);
            frac  = t_now / T_final;

            % Print every 1% of total time
            coarse = 0.01 * floor(frac / 0.01);
            if coarse > last_frac && coarse >= 0
                fprintf('Progress: %3.0f%%\n', 100 * coarse);
                last_frac = coarse;
            end

            status = 0;
            return

        otherwise
            % 'done' or other flags
            status = 0;
            return
    end
end

%% ---------------------- SINGLE-MODE a OPERATOR -----------------------

function a = annihilation_operator_sparse(nmax_local)
%ANNIHILATION_OPERATOR_SPARSE Sparse annihilation operator on Fock space 0..nmax

    vals = sqrt(1:nmax_local).';      % sqrt(n) factors
    B    = zeros(nmax_local+1, 1);
    B(1:nmax_local) = vals;

    % a|n> = sqrt(n) |n-1> -> superdiagonal in matrix (shift to lower n)
    a = spdiags(B, 1, nmax_local+1, nmax_local+1);
end

%% ---------------------- MULTI-MODE a_k OPERATORS ---------------------

function A_ops_local = multimode_annihilation_ops_sparse(nmax_local, K_local, dim_tls_local)
%MULTIMODE_ANNIHILATION_OPS_SPARSE Build {a_k} on K-mode ⊗ TLS space

    dim_ph1  = nmax_local + 1;
    a_single = annihilation_operator_sparse(nmax_local);
    Id_tls   = speye(dim_tls_local);

    A_ops_local = cell(K_local, 1);

    for kk = 1:K_local
        % Build photonic operator for mode kk as
        % I ⊗ ... ⊗ a_single (at kk) ⊗ ... ⊗ I
        op_ph = 1;
        for mm = 1:K_local
            if mm == kk
                M = a_single;
            else
                M = speye(dim_ph1);
            end
            op_ph = kron(op_ph, M);
        end

        % Extend with TLS identity
        A_ops_local{kk} = kron(op_ph, Id_tls);
    end
end

%% ---------------------- SCHRÖDINGER RHS FUNCTION ---------------------

function dpsi = schrodinger_rhs_precomputed(t, psi, H0_local, Splus_local, Sminus_local, f_handle)
%SCHRODINGER_RHS_PRECOMPUTED Right-hand side for time-dependent Schrödinger equation

    ft = f_handle(t);

    % dψ/dt = -i [ H0 + f(t) S^+ + conj(f(t)) S^- ] ψ
    dpsi = -1i * ( H0_local * psi ...
                 + ft   * (Splus_local  * psi) ...
                 + conj(ft) * (Sminus_local * psi) );
end
