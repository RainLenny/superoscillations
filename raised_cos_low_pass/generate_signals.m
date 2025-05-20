clear all

syms t w

angular_freqs = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65];
input_amps = [17.551969178838473, -49.93615200929133, 48.91892825922387, -7.222327926051708, -17.23049063370233, 8.876842152413422];

superoscillations_signal = sum(input_amps .* exp(1i * angular_freqs * t));

sinc_signal = piecewise(t == 0, 1, sinc(t/pi));

%The transition band of the filter is from omega = a to b
a = 0.65;
b = 0.85;
beta = (b-a)/(b+a);
T = 2*pi/(b+a);

filter_signal = piecewise( ...
    t == 0,                        1/T, ...
    abs(t) == T/(2*beta),          (beta*sin(pi/(2*beta)))/(2*T), ...
    sinc(t/T).*cos(pi*beta.*t/T)./(1 - (2*beta.*t/T).^2) * (1/T) ...
);
% filter_signal = piecewise(t == 0, 1, sinc(t*0.7/pi))*0.7/pi;



