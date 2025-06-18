syms t w

scaling = 1/2;

angular_freqs = [0,1,2,3]/6  * scaling;
input_amps = [-204.67096234866244, 327.0263462347299, -159.07599658671475, 37.63117153871902];

superoscillations_signal = sum(input_amps .* exp(1i * angular_freqs * t));

signal = cos(t*scaling);



