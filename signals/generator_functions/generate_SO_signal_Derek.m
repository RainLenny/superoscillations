clc;
clear;

% Add all project subfolders to search path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')));

% MATLAB Code to generate the specific superoscillating signal from the paper:
% "Reconstructing superoscillations buried deeply in noise" (Derek et al.)

% 1. Define the frequencies
% The user specified a superoscillating frequency of 1, and the constituent frequencies:
N = 4;
angular_freqs_SO = [0.3, 0.4, 0.5, 0.6];

% 2. Define the amplitudes
% Based on the specific requirements, these are the amplitudes corresponding to the frequencies above.
amps_SO = [-0.31478314809892577, 0.93822642984254, -1.0, 0.3913961476910592];

%% Save

amps_SO = amps_SO(:).';
angular_freqs_SO = angular_freqs_SO(:)';

% Save directly to the signals/data directory
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(fullfile(output_dir, 'SO_Derek.mat'), 'amps_SO', 'angular_freqs_SO');

%% PLOT

% Reconstruct the superoscillating signal over a continuous time range
t_continuous = linspace(-55, 55, 2000);
s_t = zeros(size(t_continuous));

% Reconstruct the signal: s(t) = Re( sum( c_n * exp(i * omega_n * t) ) )
for n = 1:N
    s_t = s_t + real(amps_SO(n) * exp(1i * angular_freqs_SO(n) * t_continuous));
end

% Plot the results
figure('Position', [100, 100, 800, 400]);

plot(t_continuous, s_t, 'b-', 'LineWidth', 1.5);
title('Superoscillating Function s(t) (Derek et al.)');
xlabel('Time');
xlim([-55, 55]);
grid on;

% Display the given coefficients in the command window
disp('The calculated complex coefficients c_n are:');
for i=1:N
    fprintf('c_%d = %7.3f %+.3fi\n', i, real(amps_SO(i)), imag(amps_SO(i)));
end

%% Instantaneous Frequency Analysis
% Compute instantaneous frequency if the function is available
if exist('compute_instantaneous_frequency', 'file')
    inst_freq = compute_instantaneous_frequency(s_t, t_continuous);

    figure('Color', 'w', 'Name', 'Instantaneous Frequency Analysis');
    plot(t_continuous, s_t, 'LineWidth', 2, 'DisplayName', 'wave s(t)'); 
    hold on;
    plot(t_continuous, inst_freq, 'LineWidth', 2, 'DisplayName', '\omega_{inst}(t)');
    title('SO Derek: Wave and Instantaneous Frequency');
    xlabel('Time');
    ylim([-6 6]); 
    grid on;
    legend('Location', 'northeast');

    if exist('PlotUtils', 'class')
        PlotUtils.styleAxes(gca);
    end
end
