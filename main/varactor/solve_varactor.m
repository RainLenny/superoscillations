function [tgrid, q, I, V_c] = solve_varactor(L, R, C0, alpha, V_in, tspan, y0)
% SOLVE_RLC solves the non-linear RLC circuit with varactor
% Equation: L * d^2q/dt^2 + R * dq/dt + V_c(q) = V_in(t)
% where V_c(q) = q / C0 + alpha * q^2

if nargin < 7
    y0 = [0; 0]; % [q0; I0]
end
if nargin < 4
    alpha = 0; % linear by default
end

% V_c function (Second-order nonlinearity)
Vc_func = @(q) q ./ C0 + alpha * q.^2;

% State vector y = [q; dq/dt]
% dy/dt = [y(2); (V_in(t) - R*y(2) - Vc_func(y(1))) / L]
odefun = @(t, y) [
    y(2);
    (V_in(t) - R * y(2) - Vc_func(y(1))) / L
];

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[tgrid, y] = ode45(odefun, tspan, y0, opts);

q = y(:, 1);
I = y(:, 2);
V_c = Vc_func(q);

end
