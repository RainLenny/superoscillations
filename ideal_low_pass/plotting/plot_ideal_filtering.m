generate_signals

filtered_sinc = filter_signal*pi;

% Define time spans
t_span = [-15, 15];


%% PLOT
% Create tiled layout with shared axis labels
figure;
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Subplot 1: Full range plot (t_span1)
ax1 = nexttile;
fplot(real(superoscillations_signal), t_span, 'r', 'LineWidth', 6, 'DisplayName', 'Original');
hold on;
fplot(real(superoscillations_signal), t_span, ':b', 'LineWidth', 6, 'DisplayName', 'Filtered');
% fplot(filter_signal, t_span2, 'k', 'LineWidth', 4, 'DisplayName', 'Filter');
title('Superoscillating Signal');
legend('show', 'fontweight', 'bold', 'fontsize', 12);
grid on;
xlim([-15, 15]);
ax1.FontWeight = 'bold';
ax1.FontSize = 12;
ylim([-0.5, 1]);

% Subplot 2: Zoomed-in plot (t_span2)
ax2 = nexttile;
fplot(sinc_signal, t_span, 'r', 'LineWidth', 6, 'DisplayName', 'Original');
hold on;
fplot(filtered_sinc, t_span, ':b', 'LineWidth', 6, 'DisplayName', 'Filtered');
% fplot(filter_signal, t_span1, 'k', 'LineWidth', 4, 'DisplayName', 'Filter');
title('Sinc Signal');
legend('show', 'fontweight', 'bold', 'fontsize', 12);
grid on;
xlim([-15, 15]);
ax2.FontWeight = 'bold';
ax2.FontSize = 12;
ylim([-0.5, 1]);

% Shared axis labels
xlabel(t, ['Time [2' char(960) '/' char(969) ']'], 'FontWeight', 'bold', 'FontSize', 12);
ylabel(t, 'Amplitude [arb. units]', 'FontWeight', 'bold', 'FontSize', 12);

hold off;
