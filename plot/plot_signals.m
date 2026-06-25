function [] = plot_signals(sig_configs)
% PLOT_SIGNALS Plots the time-domain and frequency-domain representations 
% of an array of signal configurations.
%
% Usage:
%   plot_signals(sig_configs)

% Add project root and all subfolders to search path
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
PlotUtils.setupDefaults();

%% Collect all frequencies to compute the fundamental period
all_freqs = [];
for i = 1:length(sig_configs)
    if isfield(sig_configs(i), 'freqs') && ~isempty(sig_configs(i).freqs)
        all_freqs = [all_freqs, sig_configs(i).freqs];
    end
end
all_freqs = unique(all_freqs);

%% SAMPLING Constants
f_sampling = 60/(2*pi);
% fundamental_period of all signals combined:
T_period = compute_fundamental_period(all_freqs, f_sampling);
signals_duration = T_period * 100; % total duration to simulate
dt = 1 / f_sampling; % sample-interval in seconds
t_axis = -signals_duration/2 : dt : signals_duration/2;
t_axis = t_axis(1:end-1);
t_axis = t_axis(:);

%% Sample the signals
for i = 1:length(sig_configs)
    sig_configs(i).sampled = sig_configs(i).data(t_axis);
end

%% Plot Time Domain
figure
hold on;
h_arr = gobjects(1, length(sig_configs));
for i = 1:length(sig_configs)
    h_arr(i) = plot(t_axis, real(sig_configs(i).sampled), '-', 'color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end
legend(h_arr);

xlabel('\boldmath$\mathbf{Time \ [2\pi/\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');

PlotUtils.styleAxes(gca);
hold off;


%% FFT the signals
N = length(t_axis);
freq_axis = linspace(-f_sampling/2, f_sampling/2, N)*2*pi; % Frequency axis

for i = 1:length(sig_configs)
    sig_configs(i).fft_mag = fftshift(abs(fft(sig_configs(i).sampled, N))) / N;
end

%% Plot FFTs
figure;
hold on;

h_arr_fft = gobjects(1, length(sig_configs));
for i = 1:length(sig_configs)
    h_arr_fft(i) = plot(freq_axis, sig_configs(i).fft_mag, '-', 'Color', sig_configs(i).color, 'DisplayName', sig_configs(i).name);
end
legend(h_arr_fft);

xlabel('\boldmath$\mathbf{Angular \ frequency \ [\omega_0]}$');
ylabel('\boldmath$\mathbf{Amplitude \ [arb]}$');

PlotUtils.styleAxes(gca);
hold off;

end