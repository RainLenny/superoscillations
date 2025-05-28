function y = sym_conv(sym_signal1, sym_signal2, t, varargin)
%SYM_CONV  Symbolic convolution of two (possibly complex) signals.
%   y = SYM_CONV(h, x, t) returns 
%       y(t) = ∫_{-∞}^{∞} h(τ) * x(t - τ) dτ
%   for symbolic h(t) and x(t).
%
%   y = SYM_CONV(h, x, t, 'causal') returns
%       y(t) = ∫_{0}^{t} h(τ) * x(t - τ) dτ
%   which is appropriate when both signals are causal (zero for t<0).
%
% Inputs:
%   sym_signal1       – a symbolic expression in t, your impulse response
%   sym_signal2       – a symbolic expression in t, your input signal
%   t       – the symbolic time variable 
%   'causal' (optional) – flag to integrate from 0 to t instead of ±∞
%
% Output:
%   y       – the symbolic result of the convolution, simplified
%
% Example:
%   syms t
%   sym_signal1 = exp(-2*t)*heaviside(t);     % causal decaying exponential
%   sym_signal2 = (sin(3*t) + 1j*cos(2*t))*heaviside(t);  % causal complex input
%   y_full   = sym_conv(h, x, t);           % full integral
%   y_causal = sym_conv(h, x, t, 'causal'); % 0→t integral
%   pretty(simplify(y_causal))
%

    % introduce dummy integration variable
    syms tau

    % substitute τ into h and (t–τ) into x
    h_tau   = subs(sym_signal1, t, tau);
    x_shift = subs(sym_signal2, t, t - tau);

    % choose limits
    if ~isempty(varargin) && strcmpi(varargin{1}, 'causal')
        a = 0; 
        b = t;
    else
        a = -inf; 
        b =  inf;
    end

    % compute and simplify
    y = simplify( int(h_tau .* x_shift, tau, a, b) );
end