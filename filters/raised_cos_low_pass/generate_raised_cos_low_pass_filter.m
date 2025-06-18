syms t

%The transition band of the filter is from omega = a to b
a = 0.7;
b = 0.8;
beta = (b-a)/(b+a);
T = 2*pi/(b+a);

filter_signal = piecewise( ...
    t == 0,                        1/T, ...
    abs(t) == T/(2*beta),          (beta*sin(pi/(2*beta)))/(2*T), ...
    sinc(t/T).*cos(pi*beta.*t/T)./(1 - (2*beta.*t/T).^2) * (1/T) ...
);




