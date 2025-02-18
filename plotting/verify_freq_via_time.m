syms t w

% Define parameters
angular_freqs = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65];
scaling_factor = 1;
angular_freqs = angular_freqs * scaling_factor;
input_amps = [17.551969178838473, -49.93615200929133, 48.91892825922387, -7.222327926051708, -17.23049063370233, 8.876842152413422];

% Define signals
superoscillations_signal = sum(input_amps .* exp(1i * angular_freqs * t));
sinc_signal = sinc(t * scaling_factor / pi);
filter_signal = sinc(t * 0.7 * scaling_factor / pi) * 0.7 / pi;

% Compute Fourier transforms
FT_superoscillations = fourier(superoscillations_signal, t, w);
FT_sinc = fourier(sinc_signal, t, w);
FT_filter = fourier(filter_signal, t, w);

% Convert symbolic expressions to functions for plotting
FT_superoscillations_fun = matlabFunction(FT_superoscillations, 'Vars', w);
FT_sinc_fun = matlabFunction(FT_sinc, 'Vars', w);
FT_filter_fun = matlabFunction(FT_filter, 'Vars', w);

% Frequency range for plotting
w_range = [-5, 5];

% Plot Fourier Transforms
figure;
subplot(3,1,1);
fplot(@(w) abs(FT_superoscillations_fun(w)), w_range);
title('Magnitude of Fourier Transform - Superoscillations Signal');
xlabel('\omega');
ylabel('|F(\omega)|');
grid on;

subplot(3,1,2);
fplot(@(w) abs(FT_sinc_fun(w)), w_range);
title('Magnitude of Fourier Transform - Sinc Signal');
xlabel('\omega');
ylabel('|F(\omega)|');
grid on;

subplot(3,1,3);
fplot(@(w) abs(FT_filter_fun(w)), w_range);
title('Magnitude of Fourier Transform - Filter Signal');
xlabel('\omega');
ylabel('|F(\omega)|');
grid on;

% Adjust layout
sgtitle('Fourier Transforms of the Signals');
