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
            if nargin < 2
                scale_y_factor = 1.0;
            end
            
            % Ensure basic properties are set correctly on the axes
            ax.TickLabelInterpreter = 'latex';
            
            % MATLAB automatically links the Legend FontSize to the Axes FontSize,
            % ignoring DefaultLegendFontSize. We manually enforce it here.
            try
                defaultLegSize = get(groot, 'DefaultLegendFontSize');
                lgds = findobj(ax.Parent, 'Type', 'legend');
                for k = 1:length(lgds)
                    % Only apply to the legend linked to this axes
                    if isequal(lgds(k).Axes, ax)
                        lgds(k).FontSize = defaultLegSize;
                    end
                end
            catch
                % If DefaultLegendFontSize is not set, do nothing
            end
            
            % Format dynamically generated ticks in bold LaTeX.
            % Double backslash (\\) is required so it isn't parsed as an escape character.
            xtickformat(ax, '$\\mathbf{%g}$');
            ytickformat(ax, '$\\mathbf{%g}$');
            
            % If y-axis is scaled, automatically place the exponent label at the top-left
            if scale_y_factor ~= 1.0
                % Disable MATLAB's native auto-exponent to prevent it from clashing 
                % visually with the custom LaTeX exponent added below.
                ax.YAxis.Exponent = 0; 
                
                k = round(log10(scale_y_factor));
                if k ~= 0
                    exponent_str = sprintf('\\times 10^{-%d}', k);
                    if k < 0
                        exponent_str = sprintf('\\times 10^{%d}', -k);
                    end
                    PlotUtils.addExponent(ax, exponent_str);
                end
            end
        end
        
        function addExponent(ax, exponent_str)
            % ADDEXPONENT Places a custom LaTeX exponent label at the top-left of the plot area.
            %   Uses normalized coordinates so (0,1) is exactly the top-left of the plot box.
            
            if nargin < 1 || isempty(ax)
                ax = gca;
            end
            
            text(ax, 0, 1, sprintf('$\\mathbf{%s}$', exponent_str), ...
                'Units', 'normalized', ...
                'Interpreter', 'latex', ...
                'FontSize', 12, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom');
        end
    end
end