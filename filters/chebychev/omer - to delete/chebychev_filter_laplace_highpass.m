function H = chebychev_filter_laplace_highpass(omega0, N, epsilon)
% computes a Chebyshev filter defined by
%
%   H(s) = prod_{k=0}^{N-1} [ v_k*(s - u_k) ] / [ u_k*(s - v_k) ],
%   places s = i*omega
%
% where
%
%   v_k = omega0 / s_k,
%   u_k = (1j*omega0) / cos( (2k+1)*pi/(2*N) ),
%
% and
%
%   s_k = -omega0*sinh((1/N)*asinh(1/epsilon))*sin((2k+1)*pi/(2*N)) ...
%         + 1j*omega0*cosh((1/N)*asinh(1/epsilon))*cos((2k+1)*pi/(2*N)).
%
% Inputs:
%   omega0  - cutoff frequency
%   N       - filter order (positive integer)
%   epsilon - ripple factor (positive scalar)
%
% Output:
%   H - the symbolic transfer function H(s)
%
% Example:
%   H = chebychev_filter(1000, 4, 0.5);
%   disp(H);

% Define the symbolic variable for s
syms s;
% Substitution of variables for a high-pass filter 
s_hp = (omega0)^2 /s;

% Precompute the constant used in the hyperbolic functions
alpha = asinh(1/epsilon);

% Initialize the filter expression
H_expr = 1;

% Loop over each stage k = 0,1,...,N-1
for k = 0:(N-1)
    % Compute theta for the current stage
    theta = (2*k + 1)*pi/(2*N);

    % Compute s_k as defined
    s_k = -omega0 * sinh(alpha/N)*sin(theta) ...
        + 1j*omega0*cosh(alpha/N)*cos(theta);

    % Compute v_k and u_k for the current stage
    v_k = omega0/s_k;
    u_k = 1j*omega0/cos(theta);

    % Multiply the filter by the current stage factor
    H_expr = H_expr * (v_k*(s_hp - u_k))/(u_k*(s_hp - v_k));
end

% Simplify the expression before returning
H = simplify(H_expr);


end
