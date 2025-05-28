clear all

syms t w

scaling = 1;

angular_freqs = [0,1,2,3]/6  * scaling;
input_amps = [-204.67096234866244, 327.0263462347299, -159.07599658671475, 37.63117153871902];

superoscillations_signal = sum(input_amps .* exp(1i * angular_freqs * t));
% superoscillations_signal =piecewise(t==0,1, real(superoscillations_signal * (sin(0.3*scaling*t)/(0.3*scaling*t))^10));

cos_signal = cos(t*scaling);

% Filter Parameters
omega0  = 1;          % Corner frequency (rad/s)
filter_order       = 7;          % Filter order
epsilon = 0.04475;    % Ripple factor


