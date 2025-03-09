generate_signals
t_span = [-15,15];

%% Plot
fig = figure;
defaultPos = get(groot, 'DefaultFigurePosition'); % Get MATLAB's default figure size

% Modify height to be half of the default
newHeight = defaultPos(4) *0.7;

% Set new figure position with the same width but half the height
set(fig, 'Position', [defaultPos(1), defaultPos(2) + defaultPos(4) / 2, defaultPos(3), newHeight]);

hold on;

fplot(real(superoscillations_signal),t_span, 'r', 'LineWidth', 6, 'DisplayName', 'Superoscilation');
fplot(sinc_signal,t_span, ':b', 'LineWidth', 6, 'DisplayName', 'Sinc');
fplot(filter_signal,t_span, '-.k', 'LineWidth', 6, 'DisplayName', 'Filter');

% Customize plot
xlabel(['Time [2' char(960) '/' char(969) ']'])
ylabel('Amplitude [arb. units]');
legend('show','fontweight','bold','fontsize',12);
grid on;
ax=gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
ylim([-0.5,1])
hold off;
