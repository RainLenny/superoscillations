clear all

filter_signal = piecewise(t == 0, 1, sinc(t*0.7/pi))*0.7/pi;



