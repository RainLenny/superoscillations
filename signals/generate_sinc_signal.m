syms t

scaling = 1;

angular_freqs = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65] * scaling;
input_amps = [17.551969178838473, -49.93615200929133, 48.91892825922387, -7.222327926051708, -17.23049063370233, 8.876842152413422];

superoscillations_signal = sum(input_amps .* exp(1i * angular_freqs * t));

signal = piecewise(t == 0, 1, sinc(t*scaling/pi));