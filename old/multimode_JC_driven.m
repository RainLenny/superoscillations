function [tgrid, Pe] = multimode_JC_driven(omega, nu0, nmax, J0, Drive_integral, T_final)
%MULTIMODE_JC_DRIVEN_IP Time-evolve a driven multi-mode JC model in the interaction picture
%
%   [tgrid, Pe] = multimode_JC_driven_IP(omega, nu0, nmax, J0, Drive_integral, T_final)
%
%   Interaction-picture dynamics with respect to the free Hamiltonian
%   H_free = sum_k ħ ω_k a_k^† a_k + (1/2)ħ ν0 σ3,
%   giving
%     V_I(t) = ħ sum_k J0 [ e^{iΔ_k t} a_k σ_+ + e^{-iΔ_k t} a_k^† σ_- ]
%              + ħ f(t) σ_+ + ħ f(t)^* σ_-,
%   where Δ_k = ν0 - ω_k and f(t) is the effective classical drive.
%
%   INPUTS:
%     omega           : 1×K vector of cavity/mode frequencies ω_k
%     nu0             : TLS transition frequency (scalar)
%     nmax            : maximum photon number per mode (0..nmax)
%     J0              : coupling strength (scalar, assumed same for all modes)
%     Drive_integral  : function handle representing the integral part in f(t),
%                       used as f(t) = -1i * J0 * Drive_integral(t) * exp(1i*nu0*t)
%                       (constants such as 𝓙_0, spectral sum, etc. can be absorbed)
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
    % omega: 1×K vector; nu0: scalar
    omega = omega(:).';                 % ensure row vector
    K     = numel(omega);               % number of modes
    Delta = nu0 - omega;                % Δ_k = ν0 - ω_k

    %% ------------------------- DIMENSIONS ------------------------------
    dim_ph1 = nmax + 1;       % local dimension per mode
    dim_ph  = dim_ph1^K;      % total photonic Hilbert space dimension
    dim_tls = 2;              % TLS dimension
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
    % A_ops{k} acts as a_k ⊗ identity on others ⊗ TLS
    A_ops = multimode_annihilation_ops_sparse(nmax, K, dim_tls);

    %% -------------- PRECOMPUTE JC MATRIX PRODUCTS (SPARSE) -------------
    % We want, for each k:
    %   J0 * e^{iΔ_k t} a_k σ_+   +   J0 * e^{-iΔ_k t} a_k^† σ_-
    % so precompute matrices:
    %   JC_plus_ops{k}  = J0 * A_k * Splus
    %   JC_minus_ops{k} = J0 * A_k^† * Sminus
    JC_plus_ops  = cell(K,1);
    JC_minus_ops = cell(K,1);
    for k = 1:K
        Ak  = A_ops{k};
        Akd = Ak';
        JC_plus_ops{k}  = J0 * (Ak  * Splus);
        JC_minus_ops{k} = J0 * (Akd * Sminus);
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
    % TLS basis assumed: [|g>; |e>] per photonic configuration
    % so excited components are entries 2,4,6,... in the full vector.
    idx_e = 2:2:dim_tot;

    %% ---------------------- DRIVE FUNCTION f(t) ------------------------
    % From your theory (with constants absorbed into Drive_integral):
    %   f(t) ∝ -i J0 ∫ dτ D^*(τ) ... e^{iν0 t}
    % We keep your earlier parametrization:
    f = @(t) -1i * J0 * Drive_integral(t) .* exp(1i * nu0 * t);

    %% --------------------- SCHRÖDINGER RHS HANDLE ---------------------
    % dψ/dt = -i V_I(t) ψ
    %
    % V_I(t) = sum_k [ J0 e^{iΔ_k t} a_k S^+ + J0 e^{-iΔ_k t} a_k^† S^- ]
    %          + f(t) S^+ + f(t)^* S^-.
    ode_rhs = @(t, psi) schrodinger_rhs_IP( ...
                            t, psi, Delta, JC_plus_ops, JC_minus_ops, ...
                            Splus, Sminus, f);

    %% ---------------------- SPARSITY PATTERN --------------------------
    % Use union of nonzero patterns as a JPattern hint
    Jpattern = spones(Splus + Sminus);
    for k = 1:K
        Jpattern = spones(Jpattern + JC_plus_ops{k} + JC_minus_ops{k});
    end

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

    % a|n> = sqrt(n) |n-1> -> superdiagonal in matrix (shift down in n)
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
        % Build photonic operator for mode kk:
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

%% --------------- SCHRÖDINGER RHS (INTERACTION PICTURE) ---------------

function dpsi = schrodinger_rhs_IP(t, psi, Delta, JC_plus_ops, JC_minus_ops, ...
                                   Splus_local, Sminus_local, f_handle)
%SCHRODINGER_RHS_IP Right-hand side for interaction-picture Schrödinger eq.
%
%   dψ/dt = -i V_I(t) ψ
%   V_I(t) = sum_k [ J0 e^{iΔ_k t} a_k S^+ + J0 e^{-iΔ_k t} a_k^† S^- ]
%            + f(t) S^+ + f(t)^* S^-.

    K  = numel(Delta);
    ft = f_handle(t);

    % Start with classical drive part
    dpsi = -1i * ( ft      * (Splus_local  * psi) ...
                 + conj(ft) * (Sminus_local * psi) );

    % Add multimode JC contributions with explicit time phases
    for k = 1:K
        dpsi = dpsi ...
             - 1i * ( exp(1i * Delta(k) * t)   * (JC_plus_ops{k}  * psi) ...
                    + exp(-1i * Delta(k) * t)  * (JC_minus_ops{k} * psi) );
    end
end
