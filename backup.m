function [tgrid, Pe] = multimode_JC_driven(omega, nu0, nmax, J_fluc, J_drive, Drive_integral, T_final)
%  Interaction-picture dynamics with respect to the free Hamiltonian.
%  Refactored for efficiency: Uses ode45 and removes unused variables.

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

% Generate Annihilation operators
A_ops  = multimode_annihilation_ops_global(basis_states, dim_tls);

%% INTERACTION OPERATORS (JAYNES–CUMMINGS TERMS)
JC_plus_ops  = cell(K,1);
JC_minus_ops = cell(K,1);

for k = 1:K
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
psi0 = kron(vac_ph, g_tls);

%% OBSERVABLE: TLS EXCITED-STATE POPULATION
% Even indices correspond to |e> in the kron(Phi, TLS) ordering
idx_e = 2:2:dim_tot;

%% CLASSICAL DRIVE
f = @(t) -1i * J_drive*J_fluc * Drive_integral(t) .* exp(1i * nu0 * t);

%% SCHRÖDINGER EQUATION (INTERACTION PICTURE)
ode_rhs = @(t,psi) schrodinger_rhs_IP( ...
    t, psi, Delta, JC_plus_ops, JC_minus_ops, Splus, Sminus, f);

%% TIME INTEGRATION
% We use ode45 (standard for unitary dynamics) and skip JPattern calculation.
outputFcn = @(t,y,flag) local_output_fun(t,y,flag,T_final);

opts = odeset('RelTol',RelTol, ...
    'AbsTol',AbsTol, ...
    'OutputFcn',outputFcn);

[tgrid, psi_all] = ode45(ode_rhs, tspan, psi0, opts);

Pe = sum(abs(psi_all(:,idx_e)).^2,2);

end

%% ======================= LOCAL FUNCTIONS =======================

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

% --- OPERATOR BUILDER (Global Cap) ---
function A_ops = multimode_annihilation_ops_global(basis, dim_tls)
    % Constructs sparse annihilation operators.
    
    dim_ph = size(basis, 1);
    K = size(basis, 2);
    Id_tls = speye(dim_tls);
    
    A_ops = cell(K,1);
    
    for k = 1:K
        row_idx = [];
        col_idx = [];
        vals    = [];
        
        for i_col = 1:dim_ph
            state_in = basis(i_col, :);
            n_k = state_in(k);
            
            if n_k > 0
                state_out = state_in;
                state_out(k) = n_k - 1;
                
                % Check index of target state (guaranteed to exist)
                [found, i_row] = ismember(state_out, basis, 'rows');
                
                if found
                    row_idx(end+1) = i_row; %#ok<AGROW>
                    col_idx(end+1) = i_col; %#ok<AGROW>
                    vals(end+1)    = sqrt(n_k); %#ok<AGROW>
                end
            end
        end
        
        Op_ph = sparse(row_idx, col_idx, vals, dim_ph, dim_ph);
        A_ops{k} = kron(Op_ph, Id_tls);
    end
end

function dpsi = schrodinger_rhs_IP(t, psi, Delta, JC_plus_ops, JC_minus_ops, Splus, Sminus, f_handle)
    Kloc = numel(Delta);
    ft = f_handle(t);

    % Drive term
    dpsi = -1i*( ft*(Splus*psi) + conj(ft)*(Sminus*psi) );

    % Interaction terms
    for kk = 1:Kloc
        p_factor = exp(1i*Delta(kk)*t);
        m_factor = conj(p_factor);
        
        dpsi = dpsi ...
            -1i*( p_factor*(JC_plus_ops{kk}*psi) ...
            + m_factor*(JC_minus_ops{kk}*psi) );
    end
end