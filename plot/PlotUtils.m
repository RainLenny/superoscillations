classdef PlotUtils
    % PlotUtils Centralized utility class for consistent publication-quality MATLAB plots.
    %   Provides methods to configure global defaults (LaTeX interpreter, fonts, grid, colors, line widths)
    %   and format axis ticks in clean bold LaTeX while preserving zoom functionality.
    
    methods (Static)
        function setupDefaults()
            % SETUPDEFAULTS Configures global graphics root defaults (groot).
            %   Sets default interpreters to LaTeX, standardizes font sizes and weights,
            %   configures line widths, sets legend options, and enables a clean white figure background.
            
            % 1. Global Interpreters: LaTeX rendering everywhere
            set(groot, 'DefaultTextInterpreter', 'latex');
            set(groot, 'DefaultLegendInterpreter', 'latex');
            set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
            
            % 2. Font Sizes and Weights
            set(groot, 'DefaultAxesFontSize', 15);
            set(groot, 'DefaultAxesFontWeight', 'bold');
            set(groot, 'DefaultTextFontSize', 14);
            set(groot, 'DefaultTextFontWeight', 'bold');
            
            % Adjust Label and Title font sizes relative to AxesFontSize (13 * 1.08 ≈ 14)
            set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.08);
            set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.08);
            set(groot, 'DefaultAxesTitleFontWeight', 'bold');
            
            % 3. Legend Defaults
            set(groot, 'DefaultLegendLocation', 'best');
            set(groot, 'DefaultLegendFontSize', 15);
            set(groot, 'DefaultLegendFontWeight', 'bold');
            
            % 4. Line & Border Defaults
            set(groot, 'DefaultLineLineWidth', 5);
            set(groot, 'DefaultAxesLineWidth', 1.5);
            
            % 5. Aesthetic Defaults: Clean white figure background
            set(groot, 'DefaultFigureColor', 'w');
            set(groot, 'DefaultFigureWindowStyle', 'docked'); % Dock figures as tabs
            
            % 6. Grid Defaults (Default to off, individual scripts can override using 'grid on')
            set(groot, 'DefaultAxesXGrid', 'off');
            set(groot, 'DefaultAxesYGrid', 'off');
        end
        
        function styleAxes(ax, scale_y_factor)
            % STYLEAXES Formats axis ticks dynamically in clean bold LaTeX and handles scaling.
            %   Applies bold LaTeX math formatting to X-ticks and Y-ticks. Optionally
            %   scales Y-ticks with a multiplier and places an appropriate exponent label.
            
            if nargin < 1 || isempty(ax)
                ax = gca;
            end
            
            auto_scale = (nargin < 2 || isempty(scale_y_factor));
            
            % Ensure basic properties are set correctly on the axes
            ax.TickLabelInterpreter = 'latex';
            ax.Box = 'off';
            
            % MATLAB automatically links the Legend FontSize to the Axes FontSize,
            % ignoring DefaultLegendFontSize. We manually enforce it here.
            try
                defaultLegSize = get(groot, 'DefaultLegendFontSize');
                lgds = findobj(ax.Parent, 'Type', 'legend');
                for k_lgd = 1:length(lgds)
                    % Only apply to the legend linked to this axes
                    if isequal(lgds(k_lgd).Axes, ax)
                        lgds(k_lgd).FontSize = defaultLegSize;
                    end
                end
            catch
                % If DefaultLegendFontSize is not set, do nothing
            end
            
            % Format dynamically generated ticks in bold LaTeX for X-axis.
            xtickformat(ax, '$\\mathbf{%g}$');
            
            % Handle Y-axis scaling and formatting per axis (supports yyaxis)
            y_axes = ax.YAxis;
            for i = 1:length(y_axes)
                y_axes(i).TickLabelFormat = '$\\mathbf{%g}$';
                
                drawnow; % Ensure tick values are populated
                
                % --- SNAP LIMITS TO TICKS ---
                % This ensures the vertical axis line ends exactly on a tick mark 
                % (giving it the clean "small horizontal line" cap at the top).
                ticks = y_axes(i).TickValues;
                if length(ticks) >= 2
                    step = ticks(2) - ticks(1);
                    lims = y_axes(i).Limits;
                    
                    % Snap min and max to the nearest outward tick step
                    new_min = floor(lims(1) / step) * step;
                    new_max = ceil(lims(2) / step) * step;
                    
                    y_axes(i).Limits = [new_min, new_max];
                    y_axes(i).TickValues = new_min:step:new_max;
                end
                % ----------------------------
                
                % Determine exponent k
                if auto_scale
                    max_val = max(abs(y_axes(i).TickValues));
                    if max_val > 0 && (max_val <= 0.1 || max_val >= 1000)
                        k = floor(log10(max_val));
                    else
                        k = 0;
                    end
                else
                    k = round(log10(scale_y_factor));
                end
                
                if k ~= 0
                    % Disable MATLAB's native auto-exponent to prevent clashing
                    y_axes(i).Exponent = 0; 
                    
                    % Setting Exponent to 0 un-scales the numbers, so we MUST 
                    % manually scale the tick values and set them as strings
                    ticks = y_axes(i).TickValues;
                    scaled_ticks = ticks / (10^k);
                    y_axes(i).TickLabels = arrayfun(@(v) sprintf('$\\mathbf{%g}$', v), scaled_ticks, 'UniformOutput', false);
                    
                    % Determine position based on which side the axis is on
                    if length(y_axes) > 1 && i == 2
                        axis_side = 'right';
                    else
                        axis_side = 'left';
                    end
                    
                    exponent_str = sprintf('\\times 10^{%d}', k);
                    PlotUtils.addExponent(ax, exponent_str, axis_side);
                else
                    % Ensure native exponent is disabled if we aren't scaling
                    y_axes(i).Exponent = 0; 
                end
            end
        end
        
        function addExponent(ax, exponent_str, axis_side)
            % ADDEXPONENT Places a custom LaTeX exponent label at the top of the plot area.
            
            if nargin < 1 || isempty(ax)
                ax = gca;
            end
            if nargin < 3
                axis_side = 'left';
            end
            
            if strcmpi(axis_side, 'right')
                x_pos = 1;
                h_align = 'right';
            else
                x_pos = 0;
                h_align = 'left';
            end
            
            text(ax, x_pos, 1.02, sprintf('$\\mathbf{%s}$', exponent_str), ...
                'Units', 'normalized', ...
                'Interpreter', 'latex', ...
                'FontSize', ax.FontSize, ...
                'HorizontalAlignment', h_align, ...
                'VerticalAlignment', 'bottom');
        end
    end
end