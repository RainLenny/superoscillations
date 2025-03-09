generate_signals

% Define time spans
t_span1 = [-70, 70];
t_span2 = [-15, 15];

%% PLOT

% Create tiled layout with shared axis labels
figure;
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Subplot 1: Full range plot (t_span1)
ax2 = nexttile;
fplot(real(superoscillations_signal), t_span2, 'r', 'LineWidth', 6, 'DisplayName', 'Superoscillation');
hold on;
fplot(sinc_signal, t_span2, ':b', 'LineWidth', 6, 'DisplayName', 'Sinc');
% fplot(filter_signal, t_span2, 'k', 'LineWidth', 4, 'DisplayName', 'Filter');
title('Zoom In');
legend('show', 'fontweight', 'bold', 'fontsize', 12);
grid on;
ax2.FontWeight = 'bold';
ax2.FontSize = 12;

% Subplot 2: Zoomed-in plot (t_span2)
ax1 = nexttile;
fplot(real(superoscillations_signal), t_span1, 'r', 'LineWidth', 6, 'DisplayName', 'Superoscillation');
hold on;
fplot(sinc_signal, t_span1, ':b', 'LineWidth', 6, 'DisplayName', 'Sinc');
% fplot(filter_signal, t_span1, 'k', 'LineWidth', 4, 'DisplayName', 'Filter');
title('Zoom Out');
legend('show', 'fontweight', 'bold', 'fontsize', 12);
grid on;
ax1.FontWeight = 'bold';
ax1.FontSize = 12;

% Shared axis labels
xlabel(t, ['Time [2' char(960) '/' char(969) ']'], 'FontWeight', 'bold', 'FontSize', 12);
ylabel(t, 'Amplitude [arb. units]', 'FontWeight', 'bold', 'FontSize', 12);

hold off;
