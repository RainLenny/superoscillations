generate_signals

f_sampling = 60/(2*pi)*10;      % sampling rate in Hz
window_length_for_filter = inf;

% Symbolic transfer function H(s) in Laplace domain
H_sym = chebychev_filter_laplace(omega0, filter_order, epsilon);

% Impulse Response h(t) via Inverse Laplace
syms s
filter_signal = ilaplace(H_sym, s, t) * heaviside(t); 

analog_filtering_func(cos_signal,superoscillations_signal,angular_freqs,...
    filter_signal, f_sampling,window_length_for_filter)