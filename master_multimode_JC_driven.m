function [tgrid, Pe, eps_trunc] = master_multimode_JC_driven(omega, nu0, nmax, J_fluc, J_drive, Drive_integral, T_final, T1, T2, do_err_est)
% MULTIMODE_JC_DRIVEN_LINDBLAD
% Solves the Lindblad Master Equation for the driven Jaynes-Cummings model
% in the interaction picture.
%
%   drho/dt = -i/hbar [V(t), rho] + sum_i L_i rho L_i' - 0.5 {L_i' L_i, rho}
%
% INPUTS:
%   omega          : 1xK vector of cavity mode frequencies
%   nu0            : TLS transition frequency
%   nmax           : Max photon number per mode
%   J_fluc         : Coupling strength (TLS-Cavity)
%   J_drive        : Coupling strength (Cavity-Drive)
%   Drive_integral : Function handle for the drive integral f(t)
%   T_final        : Simulation duration
%   T1             : Relaxation time (Spontaneous emission)
%   T2             : Dephasing time
%   do_err_est     : (Optional) Calculate truncation error
%
% OUTPUTS:
%   tgrid          : Time vector
%   Pe             : Excited state population (Trace(rho * |e><e|))
%   eps_trunc      : Truncation error estimate

%% INPUT HANDLING
if nargin < 10
    do_err_est = true;
end

%% NUMERICAL PARAMETERS
RelTol = 1e-6; % Slightly looser tol for density matrices (O(N^2))
AbsTol = 1e-8;
tspan  = [0, T_final];

%% SYSTEM PARAMETERS
omega = omega(:).';
K     = numel(omega);
Delta = nu0 - omega;

fprintf('Multimode JC Lindblad Simulation\n');
fprintf('K = %d modes, nmax = %d\n', K, nmax);
fprintf('T1 = %.2e, T2 = %.2e\n', T1, T2);

%% HILBERT SPACE DIMENSIONS
dim_ph1 = nmax + 1;
dim_ph  = dim_ph1^K;
dim_tls = 2;
dim_tot = dim_ph * dim_tls;

fprintf('Hilbert Space dim = %d (Density Matrix size = %d)\n', dim_tot, dim_tot^2);

%% LOCAL OPERATORS (TLS)
% Basis: |g> = [1;0], |e> = [0;1] is implicitly handled by the sparse construction below
% Wait, checking previous code basis: 
% vac_ph(1)=1; g_tls=[1;0] -> index 1 is |vac, g>.
% sp creates |e><g|. It maps index 1 (|g>) to index 2 (|e>).
sp = sparse([0 0; 1 0]);    % |e><g|
sm = sp.';                  % |g><e|
sz = sp*sm - sm*sp;         % |e><e| - |g><g| (Pauli Z)

%% FULL HILBERT SPACE OPERATORS
Id_ph  = speye(dim_ph);
Splus  = kron(Id_ph, sp);
Sminus = kron(Id_ph, sm);
Sz_tot = kron(Id_ph, sz);

A_ops  = multimode_annihilation_ops_sparse(nmax, K, dim_tls);

%% INTERACTION OPERATORS (HAMILTONIAN PART)
JC_plus_ops  = cell(K,1);
JC_minus_ops = cell(K,1);
for k = 1:K
    JC_plus_ops{k}  = J_fluc * (A_ops{k}  * Splus);
    JC_minus_ops{k} = J_fluc * (A_ops{k}' * Sminus);
end

%% LINDBLAD OPERATORS (DISSIPATION)
% L1 = sqrt(1/T1) * sigma_minus
% L2 = sqrt(1/2T2 - 1/4T1) * sigma_z

gamma_1 = 1/T1;
% Check for physical validity of T2
if (1/(2*T2) - 1/(4*T1)) < 0
    warning('Unphysical T2/T1 ratio. Argument of sqrt is negative. Setting dephasing to 0.');
    gamma_phi_coeff = 0;
else
    gamma_phi_coeff = sqrt(1/(2*T2) - 1/(4*T1));
end

L1 = sqrt(gamma_1) * Sminus;
L2 = gamma_phi_coeff * Sz_tot;

% Pre-compute L'L terms for the anti-commutator part
L1dagL1 = L1' * L1;
L2dagL2 = L2' * L2;

% Combine constant parts of the dissipator for efficiency?
% We keep them separate to allow clear logic in the RHS function.

%% INITIAL STATE (Density Matrix)
vac_ph = sparse(dim_ph,1);
vac_ph(1) = 1; 
g_tls = sparse([1;0]);
psi0 = kron(vac_ph, g_tls);

rho0 = psi0 * psi0'; % Pure state |psi0><psi0|
rho0_vec = rho0(:);  % Flatten to vector for ODE solver

%% CLASSICAL DRIVE FUNCTION
f = @(t) -1i * J_drive * J_fluc * Drive_integral(t) .* exp(1i * nu0 * t);

%% ODE SOLVER SETUP
% We pass the operator matrices to the local function handle
ode_rhs = @(t, rho_vec) lindblad_rhs_IP( ...
    t, rho_vec, dim_tot, Delta, JC_plus_ops, JC_minus_ops, Splus, Sminus, f, ...
    L1, L1dagL1, L2, L2dagL2);

% NOTE: JPattern for Liouvillian is dense/complex. 
% We rely on ode15s internal Jacobian estimation or switch to ode45 if not stiff.
% Given T1/T2 scales vs Rabi scales, it might be stiff.
outputFcn = @(t,y,flag) local_output_fun(t,y,flag,T_final);

opts = odeset('RelTol',RelTol, ...
    'AbsTol',AbsTol, ...
    'OutputFcn',outputFcn);

[tgrid, rho_all_vec] = ode15s(ode_rhs, tspan, rho0_vec, opts);

%% OBSERVABLES CALCULATION
% Reconstruct rho at each step and trace over cavity
nt = length(tgrid);
Pe = zeros(nt, 1);

% Operator to measure P_excited: I_ph (tensor) |e><e|
% Since Splus*Sminus = |e><e|, we can use that.
P_exc_op = Splus * Sminus; 
% Wait, Splus*Sminus = (|e><g|)(|g><e|) = |e><e|. Correct.

fprintf('Calculating observables...\n');
for i = 1:nt
    % Extract current density matrix vector and reshape
    rho_curr = reshape(rho_all_vec(i,:).', dim_tot, dim_tot);
    
    % Calculate Expectation Value: Tr(rho * Obs)
    % Efficient trace: sum(sum(rho .* Obs.'))
    Pe(i) = real(trace(rho_curr * P_exc_op));
end

%% ERROR ESTIMATION
if do_err_est
    fprintf('Calculating truncation error (running with nmax+1)...\n');
    eps_trunc = estimate_truncation_error(omega, nu0, nmax, J_fluc, J_drive, Drive_integral, T_final, T1, T2, tgrid, Pe);
    fprintf('Numerical error estimate using nmax+1: %.3e\n', eps_trunc);
else
    eps_trunc = NaN;
end

end

%% ======================= LOCAL FUNCTIONS =======================

function drho_vec = lindblad_rhs_IP(t, rho_vec, dim, Delta, JC_plus, JC_minus, Splus, Sminus, f_handle, L1, L1dagL1, L2, L2dagL2)
    
    % 1. Reshape vector to Matrix
    rho = reshape(rho_vec, dim, dim);
    
    % 2. Calculate Hamiltonian part (V)
    Kloc = numel(Delta);
    ft = f_handle(t);
    
    % V = hbar * f(t) * S_+ + h.c.
    % (hbar is absorbed in parameters)
    V = (ft * Splus) + (conj(ft) * Sminus);
    
    % Add JC terms
    for kk = 1:Kloc
        phase_p = exp(1i*Delta(kk)*t);
        phase_m = exp(-1i*Delta(kk)*t);
        V = V + (phase_p * JC_plus{kk}) + (phase_m * JC_minus{kk});
    end
    
    % Commutator -i[V, rho]
    comm_term = -1i * (V * rho - rho * V);
    
    % 3. Calculate Lindblad Dissipators
    % D[L] = L rho L' - 0.5 {L'L, rho}
    
    % Term 1: Spontaneous Emission
    D1 = (L1 * rho * L1') - 0.5 * (L1dagL1 * rho + rho * L1dagL1);
    
    % Term 2: Dephasing (only if coefficient > 0)
    % Check if L2 is zero matrix to save time (optional optimization)
    if nnz(L2) > 0
        D2 = (L2 * rho * L2') - 0.5 * (L2dagL2 * rho + rho * L2dagL2);
    else
        D2 = 0;
    end
    
    % 4. Total derivative
    drho = comm_term + D1 + D2;
    
    % 5. Flatten back to vector
    drho_vec = drho(:);
end

function eps = estimate_truncation_error(omega, nu0, nmax, J_fluc, J_drive, Drive_integral, T_final, T1, T2, t_orig, Pe_orig)
    % Recursive call with nmax + 1, do_err_est = false
    [t_new, Pe_new, ~] = master_multimode_JC_driven(omega, nu0, nmax + 1, J_fluc, J_drive, Drive_integral, T_final, T1, T2, false);
    
    % Interpolate
    Pe_new_interp = interp1(t_new, Pe_new, t_orig, 'linear');
    
    % Max difference
    eps = max(abs(Pe_new_interp - Pe_orig));
end

function a = annihilation_operator_sparse(nmax)
    vals = sqrt(1:nmax).';
    B = zeros(nmax+1,1);
    B(1:nmax) = vals;
    a = spdiags(B,1,nmax+1,nmax+1);
end

function A_ops = multimode_annihilation_ops_sparse(nmax, K, dim_tls)
    dim_ph1 = nmax + 1;
    a_single = annihilation_operator_sparse(nmax);
    A_ops = cell(K,1);
    for kk = 1:K
        op_ph = 1;
        for mm = 1:K
            if mm == kk
                M = a_single;
            else
                M = speye(dim_ph1);
            end
            op_ph = kron(op_ph,M);
        end
        A_ops{kk} = kron(op_ph, speye(dim_tls));
    end
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
            coarse = 0.1 * floor(frac/0.1); % Update every 10%
            if coarse > last_frac
                fprintf('Progress: %3.0f%%\n', 100*coarse);
                last_frac = coarse;
            end
            status = 0;
        otherwise
            status = 0;
    end
end