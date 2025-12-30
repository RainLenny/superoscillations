function [tgrid, Pe] = master_JC_drive_only(nu0, Jtot, Drive_integral, tspan, T1, T2, Ce0)
%JC_DRIVE_ONLY_MASTER Solve Master Equation dynamics for a driven JC system.
%
%   [tgrid, Pe] = JC_DRIVE_ONLY_MASTER(nu0, Jtot, Drive_integral, tspan, T1, T2, Ce0)
%
% INPUTS
%   nu0, Jtot      : System parameters for the drive (same as original)
%   Drive_integral : Function handle for the integral convolution
%   tspan          : [t0 tf] time interval
%   T1             : Relaxation time (Dissipation)
%   T2             : Dephasing time (Coherence decay)
%   Ce0            : Initial excited-state amplitude (complex). 
%                    Used to build initial density matrix. 
%                    If omitted, Ce0 = 0 (start in ground state).
%
% OUTPUTS
%   tgrid          : time points from ODE solver
%   Pe             : Population of excited state (rho_ee)
%
% DYNAMICS
%   Solves the Master Equation derived in the previous step:
%     d(rho_ee)/dt = 2*Im[ Omega * rho_eg* ] - rho_ee/T1
%     d(rho_eg)/dt = -i*Omega*(1 - 2*rho_ee) - rho_eg/T2
%   Where Omega(t) corresponds to the effective drive f(t).

    if nargin < 7 || isempty(Ce0)
        Ce0 = 0;
    end

    % --- 1. Construct Initial Density Matrix ---
    % Assume pure state initialization based on Ce0
    Cg0 = sqrt(1 - abs(Ce0)^2); 
    
    rho_ee_0 = abs(Ce0)^2;
    rho_eg_0 = Ce0 * conj(Cg0); % Coherence term

    % --- 2. Define Effective Drive (Omega) ---
    % This matches your original f(t) definition exactly
    Omega = @(t) -1i * Jtot .* Drive_integral(t) .* exp(1i * nu0 * t);

    % --- 3. Define Master Equation ODE ---
    % State vector y maps to:
    % y(1) = rho_ee (Population, real)
    % y(2) = rho_eg (Coherence, complex)
    
    odefun = @(t, y) [ ...
        % d(rho_ee)/dt
        2 * imag( Omega(t) * conj(y(2)) ) - y(1)/T1; 
        
        % d(rho_eg)/dt
        -1i * Omega(t) * (1 - 2*y(1)) - y(2)/T2
    ];

    % --- 4. Solve ---
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    
    % Note: ode45 handles the mix of real y(1) and complex y(2) automatically
    [tgrid, Y] = ode45(odefun, tspan, [rho_ee_0; rho_eg_0], opts);

    % --- 5. Output ---
    Pe = Y(:, 1); % Extract rho_ee
    
    % Ensure Pe is strictly real (remove numerical noise if any)
    Pe = real(Pe);
end