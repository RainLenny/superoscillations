syms t
scaling = 1/2;

scaling = 1;

angular_freqs = [0,1,2,3]/6  * scaling;
input_amps = [-204.67096234866244, 327.0263462347299, -159.07599658671475, 37.63117153871902];

superoscillations_signal = sum(input_amps .* exp(1i * angular_freqs * t));
superoscillations_signal =piecewise(t==0,sum(input_amps), real(superoscillations_signal * (sin(0.3*scaling*t/(2*pi))/(0.3*scaling*t/(2*pi)))^10));

signal = cos(t*scaling);


% canvas = (sin(0.3*scaling*t/(2*pi))/(0.3*scaling*t/(2*pi)))^10;
% 
% angular_freqs = [0,1,2,3]/6  * scaling;
% input_amps = [-204.67096234866244, 327.0263462347299, -159.07599658671475, 37.63117153871902];
% 
% superoscillations_signal = sum(input_amps .* exp(1i * angular_freqs * t));
% superoscillations_signal =piecewise(t==0,1, superoscillations_signal * canvas);
% 
% cos_signal = cos(t*scaling);