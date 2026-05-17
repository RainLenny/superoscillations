%% Superspectroscopy Simulation
clear; clc; close all;

%% 1. Signal Parameters (From Original Code)
T_pulse = 100;
t0 = 200;
A_SO = 7;
A_ref1 = 25;

n = 1:5;
omega = 0.18 * n; 

% Coefficients
c_SO = [-0.156 + 0.3311i, ...
        -0.862 - 1.042i, ...
         2.341 - 0.601i, ...
        -0.502 + 2.634i, ...
        -1.820 - 1.322i];

c_rms = sqrt(mean(abs(c_SO).^2));
c_flat = c_rms * ones(1, 5);

% Functions
f_SO = @(t) A_SO * sum(real(c_SO .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);
f_ref1 = @(t) A_ref1 * cos(omega(5) * t) * exp(-(t-t0)^2 / T_pulse^2);
f_ref2 = @(t) A_SO * sum(real(c_flat .* exp(1i * omega .* t))) * exp(-(t-t0)^2 / T_pulse^2);

% Generate Time-Domain Signals
t_vec = linspace(0, 600, 5000);
dt = t_vec(2) - t_vec(1);

sig_SO = zeros(size(t_vec));
sig_ref1 = zeros(size(t_vec));
sig_ref2 = zeros(size(t_vec));

for k = 1:length(t_vec)
    sig_SO(k) = f_SO(t_vec(k));
    sig_ref1(k) = f_ref1(t_vec(k));
    sig_ref2(k) = f_ref2(t_vec(k));
end

%% 2. Create Two "Almost Identical" Filters with Random Spectral Differences
% Using fir2 to design FIR filters with randomized spectral variations 
% instead of structured, monotonically increasing Butterworth cutoffs.
filter_order = 120;           % Higher order allows for finer spectral ripples
f_grid = linspace(0, 1, 500); % Normalized frequency grid (0 to 1, where 1 is Nyquist)

% Base filter magnitude profile (Low-pass covering your signal frequencies)
m1 = double(f_grid <= 0.1); 

% 1. Generate random spectral variations
rng(42);                      % Seed for consistent behavior
spectral_noise = 0.05 * randn(size(f_grid)); 

% 2. Smooth the noise using a Gaussian window to simulate realistic resonances
smoothing_window = gausswin(19); 
smoothed_noise = conv(spectral_noise, smoothing_window / sum(smoothing_window), 'same');

% 3. Create the second material's profile by adding the random ripples
m2 = m1 + smoothed_noise;
m2 = max(m2, 0);              % Ensure physical validity (magnitudes can't be < 0)

% 4. Synthesize FIR filter coefficients (FIR filters always have a denominator of 1)
b1 = fir2(filter_order, f_grid, m1);
a1 = 1;

b2 = fir2(filter_order, f_grid, m2);
a2 = 1;

% Apply filters using filtfilt (zero-phase filtering to prevent phase shifts)
E1_SO   = filtfilt(b1, a1, sig_SO);
E1_ref1 = filtfilt(b1, a1, sig_ref1);
E1_ref2 = filtfilt(b1, a1, sig_ref2);

% Transmitted signals through Material 2
E2_SO   = filtfilt(b2, a2, sig_SO);
E2_ref1 = filtfilt(b2, a2, sig_ref1);
E2_ref2 = filtfilt(b2, a2, sig_ref2);

%% 3. Calculate Distinguishability Parameter (J)
T_obs_max = 150; 
T_obs_vec = linspace(1, T_obs_max, 300);

J_SO   = zeros(size(T_obs_vec));
J_ref1 = zeros(size(T_obs_vec));
J_ref2 = zeros(size(T_obs_vec));

for k = 1:length(T_obs_vec)
    T_obs = T_obs_vec(k);
    
    % Find indices corresponding to the window [t0 - T_obs, t0 + T_obs]
    idx = find(t_vec >= (t0 - T_obs) & t_vec <= (t0 + T_obs));
    
    % Superoscillation J
    num_SO = trapz(t_vec(idx), (E1_SO(idx) - E2_SO(idx)).^2);
    den_SO = 0.5 * trapz(t_vec(idx), (E1_SO(idx).^2 + E2_SO(idx).^2));
    J_SO(k) = num_SO / den_SO;
    
    % Ref 1 J
    num_ref1 = trapz(t_vec(idx), (E1_ref1(idx) - E2_ref1(idx)).^2);
    den_ref1 = 0.5 * trapz(t_vec(idx), (E1_ref1(idx).^2 + E2_ref1(idx).^2));
    J_ref1(k) = num_ref1 / den_ref1;
    
    % Ref 2 J
    num_ref2 = trapz(t_vec(idx), (E1_ref2(idx) - E2_ref2(idx)).^2);
    den_ref2 = 0.5 * trapz(t_vec(idx), (E1_ref2(idx).^2 + E2_ref2(idx).^2));
    J_ref2(k) = num_ref2 / den_ref2;
end

%% 4. Plot Results
figure('Name', 'Superspectroscopy: Distinguishability & Filter Analysis', 'Position', [150, 150, 1050, 500]);

% --- Subplot 1: Physical Frequency Response of the Filters ---
subplot(1, 2, 1);
[h1, w] = freqz(b1, a1, 2048);
[h2, ~] = freqz(b2, a2, 2048);

% Convert digital frequency (rad/sample) to physical frequency (rad/s)
% Physical frequency = digital_frequency / dt
w_physical = w / dt; 

% Plot filter magnitudes
plot(w_physical, abs(h1), 'b-', 'LineWidth', 2); hold on;
plot(w_physical, abs(h2), 'r--', 'LineWidth', 2);

% Overlay the discrete frequencies of your input signal (omega)
for i = 1:length(omega)
    if i == 1
        xline(omega(i), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Signal Components (\omega)');
    else
        xline(omega(i), 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

title('Filter Frequency Responses (Physical Units)');
xlabel('Angular Frequency \omega (rad/s)');
ylabel('Magnitude |H(\omega)|');
grid on;

% Dynamic Limits: Zooms in comfortably to encompass the signal's full spectral window
xlim([0, omega(end) * 1.5]); 
legend('Material 1 Filter', 'Material 2 Filter', 'Location', 'best');

% --- Subplot 2: Plot Distinguishability J(T_obs) ---
subplot(1, 2, 2);
semilogy(T_obs_vec, J_SO, 'k-', 'LineWidth', 2); hold on;
semilogy(T_obs_vec, J_ref1, 'r--', 'LineWidth', 1.5);
semilogy(T_obs_vec, J_ref2, 'g:', 'LineWidth', 1.5);

title('Distinguishability (J) vs Observation Window (T_{obs})');
xlabel('Observation Window T_{obs}');
ylabel('Distinguishability Parameter J');
legend('Superoscillating (SO)', 'Ref 1: Fastest Harmonic', 'Ref 2: Flat Spectrum', 'Location', 'best');
grid on;

% Visual marker for typical superoscillatory window scale
xline(20, 'b-', 'Super-Region', 'HandleVisibility', 'off');