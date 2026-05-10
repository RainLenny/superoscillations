function [tgrid, Pe, n_k_t, n_tot_t, numerical_plus_1_error] = multimode_JC_driven_global_photon_cap(omega, nu0, nmax, J_fluc, J_drive, Drive_integral, T_final, do_err_est)
%  Interaction-picture dynamics with respect to the free Hamiltonian.
%  OPTIMIZED VERSION: Vectorized operator building + Dense state vector.

%% NUMERICAL PARAMETERS
RelTol = 1e-9;
AbsTol = 1e-9;
tspan  = [0, T_final];

%% SYSTEM PARAMETERS
omega = omega(:).';
K     = numel(omega);
Delta = nu0 - omega;

fprintf('Multimode JC simulation (Global Photon Cap)\n');
fprintf('K = %d modes, nmax (total) = %d\n', K, nmax);

%% HILBERT SPACE DIMENSIONS & BASIS
% Generate basis states where sum(n) <= nmax
basis_states = generate_fock_basis_sum_cap(nmax, K);
dim_ph = size(basis_states, 1);
dim_tls = 2;
dim_tot = dim_ph * dim_tls;

fprintf('dim_ph = %d, dim_tot = %d\n', dim_ph, dim_tot);

%% LOCAL OPERATORS (TLS)
sp = sparse([0 0; 1 0]);    % |e><g|
sm = sp.';                  % |g><e|

%% FULL HILBERT SPACE OPERATORS
Id_ph  = speye(dim_ph);

Splus  = kron(Id_ph, sp);
Sminus = kron(Id_ph, sm);

% Generate Annihilation operators (Optimized via Hashing)
A_ops  = multimode_annihilation_ops_global_fast(basis_states, dim_tls);

%% INTERACTION OPERATORS (JAYNES–CUMMINGS TERMS)
JC_plus_ops  = cell(K,1);
JC_minus_ops = cell(K,1);

for k = 1:K
    % Pre-calculate these terms to save mults inside ODE. They remain sparse.
    JC_plus_ops{k}  = J_fluc * (A_ops{k}  * Splus);
    JC_minus_ops{k} = J_fluc * (A_ops{k}' * Sminus);
end

%% INITIAL STATE
% Find the vacuum state index (all modes = 0)
vac_idx = find(sum(basis_states, 2) == 0);
if isempty(vac_idx), error('Vacuum state not found in basis'); end

vac_ph = sparse(dim_ph,1);
vac_ph(vac_idx) = 1;

g_tls = sparse([1;0]);

% OPTIMIZATION: Convert initial state to FULL vector.
% Matrices stay sparse, but the state vector becomes dense during ODE anyway.
% This prevents internal sparse-to-full thrashing in ode45.
psi0 = full(kron(vac_ph, g_tls));

%% OBSERVABLE: TLS EXCITED-STATE POPULATION
% Even indices correspond to |e> in the kron(Phi, TLS) ordering
idx_e = 2:2:dim_tot;

%% CLASSICAL DRIVE
% Pre-calculate constant factor
drive_prefactor = 1i * J_drive * J_fluc;
f = @(t) drive_prefactor * Drive_integral(t) .* exp(1i * nu0 * t);

%% SCHRÖDINGER EQUATION (INTERACTION PICTURE)
ode_rhs = @(t,psi) schrodinger_rhs_IP( ...
    t, psi, Delta, JC_plus_ops, JC_minus_ops, Splus, Sminus, f);

%% TIME INTEGRATION
outputFcn = @(t,y,flag) local_output_fun(t,y,flag,T_final);

opts = odeset('RelTol',RelTol, ...
    'AbsTol',AbsTol, ...
    'OutputFcn',outputFcn);

% ode45 is standard for non-stiff unitary dynamics
[tgrid, psi_all] = ode45(ode_rhs, tspan, psi0, opts);

Pe = sum(abs(psi_all(:,idx_e)).^2,2);

%% OBSERVABLE: PHOTON NUMBERS
% Calculate the probability of being in each photon state (tracing out TLS)
P_ph = abs(psi_all(:, 1:2:dim_tot)).^2 + abs(psi_all(:, 2:2:dim_tot)).^2;
% Expectation value of photon number for each mode <n_k>
n_k_t = P_ph * basis_states; 
% Total photon number <n_tot>
n_tot_t = sum(n_k_t, 2);

%% Numerical error calculation 
if do_err_est
    fprintf('Calculating truncation error (running with nmax+1)...\n');
    numerical_plus_1_error = estimate_truncation_error(omega, nu0, nmax, J_fluc, J_drive, Drive_integral, T_final, tgrid, Pe);
    fprintf('Numerical error estimate using nmax+1: %.3e\n', numerical_plus_1_error);
else
    % If this IS the error check run, we don't calculate an error on top of it
    numerical_plus_1_error = NaN;
end

end

%% ======================= LOCAL FUNCTIONS =======================

function eps = estimate_truncation_error(omega, nu0, nmax, J_fluc, J_drive, Drive_integral, T_final, t_orig, Pe_orig)
% The last argument 'false' prevents infinite recursion:
[t_new, Pe_new, ~, ~, ~] = multimode_JC_driven_global_photon_cap(omega, nu0, nmax + 1, J_fluc, J_drive, Drive_integral, T_final, false);
% Interpolate new result onto original time grid for comparison
Pe_new_interp = interp1(t_new, Pe_new, t_orig, 'linear');

% Compute max absolute difference
eps = max(abs(Pe_new_interp - Pe_orig));
end

function status = local_output_fun(t_curr, ~, flag, T_final)
persistent last_frac
switch flag
    case 'init'
        last_frac = -inf;
        status = 0;
    case ''
        if isempty(t_curr), status = 0; return; end
        frac = t_curr(end)/T_final;
        coarse = 0.05 * floor(frac/0.05); % Update every 5%
        if coarse > last_frac
            fprintf('Progress: %3.0f%%\n',100*coarse);
            last_frac = coarse;
        end
        status = 0;
    otherwise
        status = 0;
end
end

% --- BASIS GENERATION (Global Cap) ---
function basis = generate_fock_basis_sum_cap(nmax, K)
    % Recursively generate all vectors [n1...nK] such that sum(n) <= nmax
    if K == 1
        basis = (0:nmax)';
    else
        basis = [];
        for n = 0:nmax
            sub_basis = generate_fock_basis_sum_cap(nmax - n, K - 1);
            current_block = [n * ones(size(sub_basis, 1), 1), sub_basis];
            basis = [basis; current_block];
        end
    end
end

% --- OPTIMIZED OPERATOR BUILDER (Fast Hashing) ---
function A_ops = multimode_annihilation_ops_global_fast(basis, dim_tls)
    % Constructs sparse annihilation operators using integer hashing.
    % O(K * N * log N) instead of O(K * N^2).
    
    [dim_ph, K] = size(basis);
    Id_tls = speye(dim_tls);
    A_ops = cell(K,1);

    % Pre-compute hash weights based on the max photon number in basis
    max_n = max(basis(:));
    % Base (max_n + 1) ensures uniqueness for digits 0..max_n
    weights = (max_n + 1) .^ (0:(K-1))';
    
    % Hash the entire basis set into a single integer vector
    % basis is (dim_ph x K), weights is (K x 1) -> basis_hash is (dim_ph x 1)
    basis_hash = basis * weights;
    
    for k = 1:K
        % 1. Identify valid transitions (where n_k > 0)
        valid_mask = basis(:, k) > 0;
        src_indices = find(valid_mask);
        
        if isempty(src_indices)
            A_ops{k} = sparse(dim_ph*dim_tls, dim_ph*dim_tls);
            continue;
        end
        
        % 2. Calculate values: sqrt(n)
        % Note: The annihilation operator A acts on |n> to give sqrt(n)|n-1>
        n_vals = basis(src_indices, k);
        vals = sqrt(n_vals);
        
        % 3. Calculate target hashes
        % Decreasing photon count in mode k by 1 subtracts 1*weight(k) from hash
        target_hashes = basis_hash(src_indices) - weights(k);
        
        % 4. Find the row indices of these target hashes
        % ismember for 1D integers is drastically faster than rows
        [found, dest_indices] = ismember(target_hashes, basis_hash);
        
        % Safety check (should always find them if basis is consistent)
        if ~all(found)
            error('Error in operator construction: Target state not found.');
        end
        
        % 5. Build sparse matrix for photon space
        % Sparse takes (row, col, value)
        Op_ph = sparse(dest_indices, src_indices, vals, dim_ph, dim_ph);
        
        % 6. Tensor with TLS identity
        A_ops{k} = kron(Op_ph, Id_tls);
    end
end

function dpsi = schrodinger_rhs_IP(t, psi, Delta, JC_plus_ops, JC_minus_ops, Splus, Sminus, f_handle)
    Kloc = numel(Delta);
    ft = f_handle(t);

    % Drive term: -1i * ( f(t)*S+ + conj(f(t))*S- ) * psi
    % We compute S*psi first to avoid matrix-matrix multiplication overhead
    drive_part = ft * (Splus * psi) + conj(ft) * (Sminus * psi);

    % Interaction terms
    % Sum over k: ( exp(i*Delta*t)*JC+ + exp(-i*Delta*t)*JC- ) * psi
    
    % Accumulator for interaction part
    inter_part = zeros(size(psi));
    
    for kk = 1:Kloc
        % Calculate phase factors
        p_factor = exp(1i*Delta(kk)*t);
        
        % Apply operators to psi
        % JC_plus_ops{kk} is already (J_fluc * A * S+)
        term_plus  = p_factor * (JC_plus_ops{kk} * psi);
        term_minus = conj(p_factor) * (JC_minus_ops{kk} * psi);
        
        inter_part = inter_part + term_plus + term_minus;
    end
    
    dpsi = -1i * (drive_part + inter_part);
end