function [tgrid, Rho_ee, Rho_eg] = master_JC_drive_only(nu0, Jtot, Drive_integral, tspan, Rho0, T1, T2)
%JC_DRIVE_DISSIPATIVE Solve TLS dynamics with driving and dissipation.
%
%   [tgrid, Rho_ee, Rho_eg] = JC_DRIVE_DISSIPATIVE(nu0, Jtot, Drive_integral, tspan, Rho0, T1, T2)
%
% INPUTS
%   nu0            : TLS transition frequency
%   Jtot           : Total coupling strength (Product of J and Calligraphic J from the text)
%   Drive_integral : Function handle for the integral convolution result
%   tspan          : [t0 tf] time interval
%   Rho0           : Initial state structure (optional).
%                    Can be a scalar (Rho_ee_0) or vector [Rho_ee_0; Rho_eg_0].
%                    Default is ground state (Rho_ee=0).
%   T1             : Spontaneous emission lifetime (set to inf if neglected)
%   T2             : Coherence time (set to inf if neglected)
%
% OUTPUTS
%   tgrid          : time points from ODE solver
%   Rho_ee         : Excited state population (rho_ee) over time
%   Rho_eg         : Coherence (rho_eg) over time
%
% PHYSICS CONTEXT (from User LaTeX)
%   Solves the simplified Lindblad master equation:
%   d(rho_ee)/dt = 2 Im[ V(t) * rho_eg* ] - rho_ee/T1
%   d(rho_eg)/dt = -i V(t) * (1 - 2*rho_ee) - rho_eg/T2
%
%   Where V(t) matches the effective drive f(t) in the previous code.

% --- Default Argument Handling ---
if nargin < 5 || isempty(Rho0)
    Rho_ee_0 = 0;
    Rho_eg_0 = 0;
elseif isscalar(Rho0)
    Rho_ee_0 = Rho0;
    Rho_eg_0 = 0; % Assume no initial coherence if only population is given
else
    Rho_ee_0 = Rho0(1);
    Rho_eg_0 = Rho0(2);
end

if nargin < 6 || isempty(T1), T1 = inf; end
if nargin < 7 || isempty(T2), T2 = inf; end

% --- Effective Drive Definition ---
% Matching the variable 'f' from your previous code.
% In your LaTeX, V(t) = hbar * J * J_cal * f(t).
% Here, we assume Jtot accounts for the coupling constants.
% f_drive represents the matrix element V_eg(t)/hbar.
f_drive = @(t) -1i .* Drive_integral(t) .* exp(1i * nu0 * t);

% --- Lindblad ODE System ---
% y(1) = rho_ee (Real)
% y(2) = rho_eg (Complex)

odefun = @(t, y) [ ...
    % d(rho_ee)/dt = 2 * Im( f * conj(rho_eg) ) - rho_ee/T1
    2 *Jtot* imag( f_drive(t) * conj(y(2)) ) - y(1)/T1;

    % d(rho_eg)/dt = -i * f * (1 - 2*rho_ee) - rho_eg/T2
    -1i *Jtot* f_drive(t) * (1 - 2*y(1)) - y(2)/T2
    ];

% --- Solver Configuration ---
% Note: ode45 handles mixed real/complex states automatically.
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

[tgrid, y] = ode45(odefun, tspan, [Rho_ee_0; Rho_eg_0], opts);

% --- Output Formatting ---
Rho_ee = y(:, 1); % Population
Rho_eg = y(:, 2); % Coherence
end