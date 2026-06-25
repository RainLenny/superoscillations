function plot_instantaneous_frequency(sig_configs, t_axis, x_limits, y_limits)
% PLOT_INSTANTANEOUS_FREQUENCY Plots the wave and instantaneous frequency.
%
%   plot_instantaneous_frequency(sig_configs, t_axis)
%   plot_instantaneous_frequency(sig_configs, t_axis, x_limits)
%   plot_instantaneous_frequency(sig_configs, t_axis, x_limits, y_limits)

    if nargin < 3
        x_limits = [];
    end
    if nargin < 4
        y_limits = [];
    end

    dt = t_axis(2) - t_axis(1);
    num_sigs = length(sig_configs);
    
    figure('Color', 'w', 'Name', 'Instantaneous Frequency Analysis');

    for i = 1:num_sigs
        % Handle cases where the signal is provided as a function handle
        if isa(sig_configs(i).data, 'function_handle')
            y_sig = arrayfun(sig_configs(i).data, t_axis);
        else
            y_sig = sig_configs(i).data;
        end
        
        inst_freq = compute_instantaneous_frequency(y_sig, dt);
        
        subplot(num_sigs, 1, i);
        hold on;
        
        if isfield(sig_configs, 'color') && ~isempty(sig_configs(i).color)
            plot_color = sig_configs(i).color;
        else
            plot_color = 'b';
        end
        
        plot(t_axis, real(y_sig), 'LineWidth', 2, 'Color', plot_color, 'DisplayName', 'wave');
        plot(t_axis, inst_freq, 'LineWidth', 2, 'Color', [0.85, 0.33, 0.1], 'DisplayName', 'd\_angle/dt');
        
        if isfield(sig_configs, 'name')
            title(sprintf('%s : Wave and Instantaneous Frequency', strrep(sig_configs(i).name, '\', '\\')));
        else
            title(sprintf('Signal %d : Wave and Instantaneous Frequency', i));
        end
        
        xlabel('time');
        
        if ~isempty(x_limits)
            xlim(x_limits);
        end
        if ~isempty(y_limits)
            ylim(y_limits);
        end
        
        grid on;
        legend('Location', 'northeast');
        
        if exist('PlotUtils', 'class')
            PlotUtils.styleAxes(gca);
        end
    end
end
