function [t, rho] = lambda_system_dynamics(t_span, c_init, Omega_p_func, Omega_s_func, delta, Delta)
% LAMBDA_SYSTEM_DYNAMICS Solves the 3-level lambda system dynamics
%
% H(t) = (1/2) * [ 0,          0,           Omega_p(t);
%                  0,          -2*delta,    Omega_s(t);
%                  Omega_p(t)', Omega_s(t)', -2*Delta ]
%
% c_dot = -i H(t) c
%
% INPUTS:
%   t_span       - [t_start t_end] or time vector for ode45
%   c_init       - Initial state vector [c1_0; c2_0; c3_0]
%   Omega_p_func - Function handle for Probe Rabi frequency Omega_p(t)
%   Omega_s_func - Function handle for Stokes Rabi frequency Omega_s(t)
%   delta        - Two-photon detuning
%   Delta        - Single-photon detuning
%
% OUTPUTS:
%   t            - Time points
%   rho          - Populations [|c1|^2, |c2|^2, |c3|^2]

    dydt = @(t, c) -1i * 0.5 * [ ...
        0, 0, Omega_p_func(t); ...
        0, -2*delta, Omega_s_func(t); ...
        conj(Omega_p_func(t)), conj(Omega_s_func(t)), -2*Delta ...
    ] * c;

    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    [t, c_out] = ode45(dydt, t_span, c_init, options);

    % Populations
    rho = abs(c_out).^2;
end
