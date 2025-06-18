f_sampling = 60/(2*pi)*10;      % sampling rate in Hz
window_length_for_filter = inf;

% Filter Parameters
omega0  = 1;          % Corner frequency (rad/s)
filter_order       = 7;          % Filter order
epsilon = 0.04475;    % Ripple factor

% Symbolic transfer function H(s) in Laplace domain
H_sym = chebychev_filter_laplace(omega0, filter_order, epsilon);

% Impulse Response h(t) via Inverse Laplace
syms s
filter_signal = ilaplace(H_sym, s, t) * heaviside(t); 

