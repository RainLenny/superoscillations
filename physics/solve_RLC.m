function [tgrid, q, I, E_L, E_C, E_tot] = solve_RLC(L, R, C0, alpha, V_in, tspan, y0)
% SOLVE_RLC solves the non-linear RLC circuit with varactor
% Equation: L * d^2q/dt^2 + R * dq/dt + V_c(q) = V_in(t)
% where V_c(q) = q / C0 + alpha * q^3

if nargin < 7
    y0 = [0; 0]; % [q0; I0]
end
if nargin < 4
    alpha = 0; % linear by default
end

% V_c function
Vc = @(q) q ./ C0 + alpha * q.^3;

% Energy of the capacitor
Ec_func = @(q) (q.^2) ./ (2*C0) + (alpha * q.^4) ./ 4;

% State vector y = [q; dq/dt]
% dy/dt = [y(2); (V_in(t) - R*y(2) - Vc(y(1))) / L]
odefun = @(t, y) [
    y(2);
    (V_in(t) - R * y(2) - Vc(y(1))) / L
];

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[tgrid, y] = ode45(odefun, tspan, y0, opts);

q = y(:, 1);
I = y(:, 2);

E_L = 0.5 * L * I.^2;
E_C = Ec_func(q);
E_tot = E_L + E_C;

end
