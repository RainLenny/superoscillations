function H = chebychev_filter_magnitude(omega0, N, epsilon)
% returns a symbolic expression for the squared magnitude 
% response of a Chebyshev filter.
%
%   |H(omega)|^2 = (epsilon^2 * T_N^2(omega0/omega)) / (1 + epsilon^2 * T_N^2(omega0/omega))
%
% Inputs:
%   N       - Order of the filter (non-negative integer)
%   omega0  - Corner frequency
%   epsilon - Ripple factor
%
% Output:
%   H2_expr - A symbolic expression in the variable omega representing
%             the squared magnitude response.
%
% Example:
%   syms omega
%   H2_expr = chebyH2_sym(3, 1, 0.5);
%   pretty(H2_expr);
%

    % Define the symbolic variable
    syms omega
    
    % Compute the argument for the Chebyshev polynomial
    x = omega0/omega;
    
    % Compute the Chebyshev polynomial T_N(x) symbolically.
    Tn = chebyshevT(N, x);
    
    % Form the squared magnitude expression.
    H = sqrt((epsilon^2 * Tn^2) / (1 + epsilon^2 * Tn^2));
    % H = sqrt(1 / (1 + epsilon^2 * Tn^2));

end
