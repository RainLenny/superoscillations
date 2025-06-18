function y = sym_conv(sym_signal1, sym_signal2, t, varargin)
%SYM_CONV  Symbolic convolution via Fourier-domain multiplication.
%   y = SYM_CONV(h, x, t) returns
%       y(t) = ifourier(fourier(h(t)) * fourier(x(t)))
%   y = SYM_CONV(h, x, t, 'causal') enforces causality by zeroing
%       out y(t) for t < 0 (i.e. multiplying by heaviside(t)).
%
% Inputs:
%   sym_signal1    – symbolic expression in t (e.g., impulse response h(t))
%   sym_signal2    – symbolic expression in t (e.g., input signal x(t))
%   t              – symbolic time variable
%   'causal'       – optional flag to enforce y(t)=0 for t<0
%
% Output:
%   y              – symbolic expression for the convolution result, simplified

    % define frequency variable
    syms omega real

    % compute Fourier transforms
    H = fourier(sym_signal1, t, omega);
    X = fourier(sym_signal2, t, omega);

    % multiply in frequency domain and invert
    Y = simplify(H .* X);
    y_full = simplify(ifourier(Y, omega, t));

    % apply causality if requested
    if ~isempty(varargin) && strcmpi(varargin{1}, 'causal')
        y = simplify(y_full * heaviside(t));
    else
        y = y_full;
    end
end
