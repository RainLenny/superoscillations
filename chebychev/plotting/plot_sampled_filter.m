%% Plot Magnitude and Phase Response of a 7th Order Low-Pass Chebyshev Filter

clear; clc; close all;

%% Filter Parameters
omega0 = 1;  % Corner frequency (rad/s)
filter_order = 7;         % Filter order
epsilon = 0.04475;  % Ripple factor 

H_sym = chebychev_filter_magnitude(omega0, filter_order, epsilon);

% Convert the symbolic function to a numeric function handle.
H_func = matlabFunction(H_sym, 'Vars', sym('omega'));

%% Frequency Range and Sampling
w_min = -30;
w_max = 30;
num_samples = 50000;                     % Number of sample points
w = linspace(w_min, w_max, num_samples); % Frequency vector

%% Evaluate the Filter Response at Sampled Frequencies
H_samples = H_func(w);

%% Calculate Impulse Response via Inverse FFT
% Note: The impulse response here is computed via an IFFT of the finite frequency 
% response samples. Since the frequency response is only sampled over [w_min, w_max],
% the resulting impulse response is an approximation.

num_samples = length(H_samples)*10;

% Desired IFFT length with zero padding
filter_order = num_samples;

% Frequency spacing remains the same:
delta_w = (w_max - w_min) / (num_samples - 1);

% The total time span T (same as before):
T = 2*pi / delta_w;

% Create a time vector with N points 
t = linspace(-T*(filter_order-1)/(2*filter_order), T*(filter_order-1)/(2*filter_order), filter_order);

% Now perform the zero-padded IFFT
h_time = fftshift(ifft(ifftshift(H_samples)));


%% Plot the Frequency Response
figure;

% Plot Magnitude Response in dB
subplot(2,1,1);
plot(w, abs(H_samples), '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Magnitude Response');
hold on;
legend('show');
title('Magnitude Response of the Chebyshev Filter (dB)');
xlabel('Frequency (Rad/s)');
ylabel('Magnitude');
grid on;

% Plot Phase Response
subplot(2,1,2);
plot(w, unwrap(angle(H_samples)), '-o', 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Phase Response');
hold on;
legend('show');
title('Phase Response of the Chebyshev Filter');
xlabel('Frequency (Rad/s)');
ylabel('Phase');
grid on;

%% Plot the Impulse Response in Time Domain
figure;
plot(real(h_time), '-o', 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Impulse Response');
hold on;
legend('show');
title('Impulse Response of the Chebyshev Filter');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

