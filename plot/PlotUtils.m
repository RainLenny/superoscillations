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
        
        function styleAxes(ax, scale_y_factor, scale_x_factor)
            % STYLEAXES Formats axis ticks dynamically in clean bold LaTeX and handles scaling.
            %   Applies bold LaTeX math formatting to X-ticks and Y-ticks. Optionally
            %   scales Y-ticks and X-ticks with a multiplier and places an appropriate exponent label.
            
            if nargin < 1 || isempty(ax)
                ax = gca;
            end
            
            auto_scale_y = (nargin < 2 || isempty(scale_y_factor));
            auto_scale_x = (nargin < 3 || isempty(scale_x_factor));
            
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
            
            % Handle X-axis scaling and formatting
            x_axes = ax.XAxis;
            for i = 1:length(x_axes)
                x_axes(i).TickLabelFormat = '$\\mathbf{%g}$';
                drawnow;
                
                % --- SNAP LIMITS TO TICKS ---
                if strcmp(x_axes(i).LimitsMode, 'auto')
                    ticks = x_axes(i).TickValues;
                    if length(ticks) >= 2
                        step = ticks(2) - ticks(1);
                        lims = x_axes(i).Limits;
                        new_min = floor(lims(1) / step) * step;
                        new_max = ceil(lims(2) / step) * step;
                        x_axes(i).Limits = [new_min, new_max];
                        x_axes(i).TickValues = new_min:step:new_max;
                    end
                end
                
                if auto_scale_x
                    max_val = max(abs(x_axes(i).TickValues));
                    if max_val > 0 && (max_val <= 0.1 || max_val >= 1000)
                        k_x = floor(log10(max_val));
                    else
                        k_x = 0;
                    end
                else
                    k_x = round(log10(scale_x_factor));
                end
                
                if k_x ~= 0
                    x_axes(i).Exponent = 0;
                    ticks = x_axes(i).TickValues;
                    scaled_ticks = ticks / (10^k_x);
                    x_axes(i).TickLabels = arrayfun(@(v) sprintf('$\\mathbf{%g}$', v), scaled_ticks, 'UniformOutput', false);
                    exponent_str = sprintf('\\times 10^{%d}', k_x);
                    PlotUtils.addExponent(ax, exponent_str, 'x');
                else
                    x_axes(i).Exponent = 0;
                end
            end
            
            % Handle Y-axis scaling and formatting per axis (supports yyaxis)
            y_axes = ax.YAxis;
            for i = 1:length(y_axes)
                y_axes(i).TickLabelFormat = '$\\mathbf{%g}$';
                
                drawnow; % Ensure tick values are populated
                
                % --- SNAP LIMITS TO TICKS ---
                % This ensures the vertical axis line ends exactly on a tick mark 
                % (giving it the clean "small horizontal line" cap at the top).
                if strcmp(y_axes(i).LimitsMode, 'auto')
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
                end
                % ----------------------------
                
                % Determine exponent k
                if auto_scale_y
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
        
        function addExponent(ax, exponent_str, axis_loc)
            % ADDEXPONENT Places a custom LaTeX exponent label.
            
            if nargin < 1 || isempty(ax)
                ax = gca;
            end
            if nargin < 3
                axis_loc = 'left';
            end
            
            switch lower(axis_loc)
                case 'right'
                    try
                        x_pos = ax.YAxis(2).Label.Position(1);
                    catch
                        x_pos = ax.XLim(2);
                    end
                    y_pos = ax.YLim(2);
                    h_align = 'center';
                    v_align = 'bottom';
                    txt_units = 'data';
                case 'left'
                    x_pos = ax.YAxis(1).Label.Position(1);
                    y_pos = ax.YLim(2);
                    h_align = 'center';
                    v_align = 'bottom';
                    txt_units = 'data';
                case 'x'
                    % Position the exponent at the same horizontal line as the x-axis title
                    % using data units so it tracks dynamically when figure is resized/exported.
                    x_pos = ax.XLim(2);
                    y_pos = ax.XLabel.Position(2);
                    h_align = 'right';
                    v_align = ax.XLabel.VerticalAlignment;
                    txt_units = 'data';
                otherwise
                    x_pos = 0;
                    y_pos = 1.02;
                    h_align = 'left';
                    v_align = 'bottom';
                    txt_units = 'normalized';
            end
            
            txt = text(ax, x_pos, y_pos, sprintf('$\\mathbf{%s}$', exponent_str), ...
                'Units', txt_units, ...
                'Interpreter', 'latex', ...
                'FontSize', max(1, ax.FontSize ), ...
                'HorizontalAlignment', h_align, ...
                'VerticalAlignment', v_align);
                
            if strcmpi(axis_loc, 'x') || strcmpi(axis_loc, 'left') || strcmpi(axis_loc, 'right')
                % Add MarkedClean listener to perfectly track layout changes (like label additions or export resizes)
                addlistener(ax, 'MarkedClean', @(src, ev) PlotUtils.updateExponentPos(txt, ax, axis_loc));
            end
        end
        
        function updateExponentPos(txt, ax, axis_loc)
            % UPDATEEXPONENTPOS Safely updates the position of the exponent labels
            if isvalid(txt) && isvalid(ax)
                try
                    if strcmpi(axis_loc, 'x')
                        new_pos = [ax.XLim(2), ax.XLabel.Position(2), 0];
                    elseif strcmpi(axis_loc, 'left')
                        new_pos = [ax.YAxis(1).Label.Position(1), ax.YLim(2), 0];
                    elseif strcmpi(axis_loc, 'right')
                        try
                            new_pos = [ax.YAxis(2).Label.Position(1), ax.YLim(2), 0];
                        catch
                            new_pos = [ax.XLim(2), ax.YLim(2), 0];
                        end
                    else
                        return;
                    end
                    
                    if any(txt.Position ~= new_pos)
                        txt.Position = new_pos;
                    end
                catch
                end
            end
        end
    end
end