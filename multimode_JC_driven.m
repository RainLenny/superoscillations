function [tgrid, Pe, eps_trunc] = multimode_JC_driven(omega, nu0, nmax, J, Drive_integral, T_final, do_err_est)
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
%     J              : coupling strength (scalar, assumed same for all modes)
%     Drive_integral  : function handle representing the integral part in f(t),
%                       used as f(t) = -1i * J * Drive_integral(t) * exp(1i*nu0*t)
%     T_final         : final simulation time
%
%   OUTPUTS:
%     tgrid : time points from ODE solver
%     Pe    : excited-state population of the TLS at each time in tgrid
%   do_err_est : (Optional) Boolean. If true (default), calculates truncation
%                error by running the simulation again with nmax+1.

%% INPUT HANDLING
if nargin < 7
    do_err_est = true; % Default to calculating error
end

%% NUMERICAL PARAMETERS for the ODE solver
RelTol = 1e-9;
AbsTol = 1e-9;
tspan  = [0, T_final];

%% SYSTEM PARAMETERS
omega = omega(:).';
K     = numel(omega);
Delta = nu0 - omega;


fprintf('Multimode JC simulation\n');
fprintf('K = %d modes, nmax = %d\n', K, nmax);

%% HILBERT SPACE DIMENSIONS
dim_ph1 = nmax + 1;
dim_ph  = dim_ph1^K;
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

A_ops  = multimode_annihilation_ops_sparse(nmax, K, dim_tls);

%% INTERACTION OPERATORS (JAYNES–CUMMINGS TERMS)
JC_plus_ops  = cell(K,1);
JC_minus_ops = cell(K,1);

for k = 1:K
    JC_plus_ops{k}  = J * (A_ops{k}  * Splus);
    JC_minus_ops{k} = J * (A_ops{k}' * Sminus);
end

%% INITIAL STATE
vac_ph = sparse(dim_ph,1);
vac_ph(1) = 1;
g_tls = sparse([1;0]);
psi0 = kron(vac_ph, g_tls);

%% OBSERVABLE: TLS EXCITED-STATE POPULATION
idx_e = 2:2:dim_tot;

%% CLASSICAL DRIVE
f = @(t) -1i * J * Drive_integral(t) .* exp(1i * nu0 * t);

%% SCHRÖDINGER EQUATION (INTERACTION PICTURE)
% We pass the operator matrices to the local function handle
ode_rhs = @(t,psi) schrodinger_rhs_IP( ...
    t, psi, Delta, JC_plus_ops, JC_minus_ops, Splus, Sminus, f);

%% JACOBIAN SPARSITY PATTERN
Jpattern = spones(Splus + Sminus);
for k = 1:K
    Jpattern = spones(Jpattern + JC_plus_ops{k} + JC_minus_ops{k});
end

%% TIME INTEGRATION

outputFcn = @(t,y,flag) local_output_fun(t,y,flag,T_final);

opts = odeset('RelTol',RelTol, ...
    'AbsTol',AbsTol, ...
    'JPattern',Jpattern, ...
    'OutputFcn',outputFcn);

[tgrid, psi_all] = ode15s(ode_rhs, tspan, psi0, opts);

Pe = sum(abs(psi_all(:,idx_e)).^2,2);

%% ERROR ESTIMATION (RECURSIVE CALL)
if do_err_est
    fprintf('Calculating truncation error (running with nmax+1)...\n');
    eps_trunc = estimate_truncation_error(omega, nu0, nmax, J, Drive_integral, T_final, tgrid, Pe);
    fprintf('Numerical error estimate using nmax+1: %.3e\n', eps_trunc);
else
    % If this IS the error check run, we don't calculate an error on top of it
    eps_trunc = NaN;
end

end

%% ======================= LOCAL FUNCTIONS =======================

function eps = estimate_truncation_error(omega, nu0, nmax, J0, Drive_integral, T_final, t_orig, Pe_orig)

% The last argument 'false' prevents infinite recursion:
[t_new, Pe_new, ~] = multimode_JC_driven(omega, nu0, nmax + 1, J0, Drive_integral, T_final, false);

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
        coarse = 0.01 * floor(frac/0.01);
        if coarse > last_frac
            fprintf('Progress: %3.0f%%\n',100*coarse);
            last_frac = coarse;
        end
        status = 0;
    otherwise
        status = 0;
end
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
Id_tls = speye(dim_tls);

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
    A_ops{kk} = kron(op_ph,Id_tls);
end
end

function dpsi = schrodinger_rhs_IP(t, psi, Delta, JC_plus_ops, JC_minus_ops, Splus, Sminus, f_handle)
Kloc = numel(Delta);
ft = f_handle(t);

dpsi = -1i*( ft*(Splus*psi) + conj(ft)*(Sminus*psi) );

for kk = 1:Kloc
    dpsi = dpsi ...
        -1i*( exp(1i*Delta(kk)*t)*(JC_plus_ops{kk}*psi) ...
        + exp(-1i*Delta(kk)*t)*(JC_minus_ops{kk}*psi) );
end
end