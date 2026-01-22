%% Setup Parameters
clear; clc; close all;

% --- INPUTS ---
K = 5;                  % System parameter
t = 500;                % Time parameter
J_values = [0.001,0.0015];    % Different J values to compare
n_max_range = 0:1:500;   % X-axis: Range of n_max (truncation point)
% --------------

%% Calculation Loop
results = nan(length(n_max_range), length(J_values));

for j_idx = 1:length(J_values)
    current_J = J_values(j_idx);
    for n_idx = 1:length(n_max_range)
        n_max = n_max_range(n_idx);
        
        val = calculate_pe_bound(K, current_J, t, n_max);
        
        % Filter: If value drops below machine epsilon, set to NaN to avoid
        % plotting noise at the bottom of the log chart.
        if val < eps
            results(n_idx, j_idx) = NaN; 
        else
            results(n_idx, j_idx) = val;
        end
    end
end

%% Generate Proper Log Plot
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% 'semilogy' automatically sets the Y-axis to log scale
h = semilogy(n_max_range, results, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5);

% --- Proper Log Scale Formatting ---

% 1. Grid Minor: Essential for reading log plots
grid on;
grid minor; 

% 2. Axes Labels
xlabel('Truncation Order (n_{max})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error Bound |\Delta P_e(t)|', 'FontSize', 12, 'FontWeight', 'bold');
title({['Error Bound Convergence for K=' num2str(K) ', t=' num2str(t)]; ...
       '(Logarithmic Scale)'}, 'FontSize', 14);

% 3. Y-Axis Limits
% Lock the view to relevant data. Because n! grows fast, this curve drops
% steeply. We limit the bottom to 10^-16 (machine precision).
ylim([1e-16, max(results(:)) * 5]); 

% 4. Dynamic Legend
legend_labels = arrayfun(@(x) sprintf('J = %.3f', x), J_values, 'UniformOutput', false);
legend(h, legend_labels, 'Location', 'SouthWest', 'FontSize', 11);

% 5. Aesthetic tweaks
set(gca, 'FontSize', 11, 'LineWidth', 1.2); % Make tick marks clearer