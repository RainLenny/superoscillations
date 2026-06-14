function [t, rho] = optical_bloch(t_span, rho_init, nu0, params, f_func)
% SOLVE_BLOCH_ODES Solves:
%   rho1_dot =  nu0*rho2 - rho1/T2
%   rho2_dot = -nu0*rho1 + 2*OmegaTilde*f(t)*rho3 - rho2/T2
%   rho3_dot = -2*OmegaTilde*f(t)*rho2 - (rho3 - rho30)/T1
%
% INPUTS:
%   t_span      - [t_start t_end] or time vector for ode45
%   rho_init    - [rho1_0 rho2_0 rho3_0]
%   nu0         - nu_0
%   OmegaTilde  - \tilde{\Omega}
%   params      - struct with fields: T1, T2, rho30
%   f_func      - function handle for f(t), e.g. @(t) 1 or @(t) sin(t)
%
% OUTPUTS:
%   t           - time points
%   rho         - [rho1 rho2 rho3] columns

T1  = params.T1;
T2  = params.T2;
r30 = params.rho30;
OmegaTilde = params.OmegaTilde;

if T1 == 0 || T2 == 0
    error('T1 and T2 must be non-zero.');
end

% y(1)=rho1, y(2)=rho2, y(3)=rho3
dydt = @(t, y) [ ...
    nu0 * y(2) - y(1)/T2; ...
    -nu0 * y(1) + 2*OmegaTilde * f_func(t) * y(3) - y(2)/T2; ...
    -2*OmegaTilde * f_func(t) * y(2) - (y(3) - r30)/T1 ...
    ];

% Tighten tolerances to strictly preserve the Bloch vector length.
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[t, rho] = ode45(dydt, t_span, rho_init, options);
end