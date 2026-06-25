function sig_configs = prepare_signal_config(sig_configs)
% PREPARE_SIGNAL_CONFIG Validates a struct array of signal configurations
% and assigns distinct default colors for any signals that lack one.
%
% Usage:
%   sig_configs = prepare_signal_config(sig_configs);
%
% Inputs:
%   sig_configs - An array of structs. Each struct must have 'name' and 'data'.
%                 Optionally, it can have 'color', 'freqs', or other fields.
%
% Outputs:
%   sig_configs - The validated struct array with 'color' assigned where missing.

    % Default color palette (distinct colors)
    default_colors = {
        'r',            % Red
        'b',            % Blue
        [0, 0.5, 0],    % Dark Green
        'm',            % Magenta
        [0.85, 0.33, 0.1], % Orange
        'c',            % Cyan
        'k',            % Black
        [0.49, 0.18, 0.56], % Purple
        [0.93, 0.69, 0.13]  % Yellow/Gold
    };

    num_colors = length(default_colors);
    color_idx = 1; % Start at the first default color

    for i = 1:length(sig_configs)
        % Ensure standard fields exist if they don't already
        if ~isfield(sig_configs, 'name')
            error('Signal at index %d must have a ''name'' field.', i);
        end
        if ~isfield(sig_configs, 'data')
            error('Signal at index %d must have a ''data'' field (function handle).', i);
        end
        if ~isfield(sig_configs, 'color')
            sig_configs(i).color = '';
        end

        % Assign color if empty
        if isempty(sig_configs(i).color)
            sig_configs(i).color = default_colors{mod(color_idx - 1, num_colors) + 1};
            color_idx = color_idx + 1;
        end
    end
end
