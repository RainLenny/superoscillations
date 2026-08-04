function [tgrid, Pe, Ce, Cg] = JC_drive_only(nu0, Jtot, Drive_integral, tspan, Ce0)
%JC_DRIVE_ONLY Solve two-level dynamics for a driven JC system.
%
%   [tgrid, Pe, Ce, Cg] = JC_DRIVE_ONLY(nu0, Jtot, Drive_integral, tspan, Ce0)
%
% INPUTS
%The Drive_integral should be the integral convolution result 
%from the effective drive calculation, the other constants and nu0
%phase are had here
%   tspan  : [t0 tf] time interval
%   Ce0    : initial excited-state amplitude (complex). If omitted,
%            Ce0 = 0 (start in ground state).
%
% OUTPUTS
%   tgrid  : time points from ODE solver
%   Pe     : |Ce(t)|^2 excitation probability
%
% NOTES
%   Uses effective drive  f(t) = i J0 Drive_integral exp(i nu0 t),
%   where J0 and nu0 are assumed to be defined in the workspace.

    if nargin < 5 || isempty(Ce0)
        Ce0 = 0;
    end

    Cg0 = sqrt(1 - abs(Ce0)^2);

    % Effective drive
    f = @(t) -1i * Drive_integral(t) .* exp(1i * nu0 * t);

    % Two-level ODE: d/dt [Ce; Cg] = -i [0 f; f* 0] [Ce; Cg]
    odefun = @(t, C) [ ...
        -1i * Jtot * f(t)        * C(2);  % dCe/dt
        -1i * Jtot * conj(f(t))  * C(1)   % dCg/dt
    ];

    opts = odeset('RelTol', 100 * eps, 'AbsTol', 1e-20);
    [tgrid, C] = ode45(odefun, tspan, [Ce0; Cg0], opts);

    Pe = abs(C(:, 1)).^2;
    Ce = C(:, 1);
    Cg = C(:, 2);
end
