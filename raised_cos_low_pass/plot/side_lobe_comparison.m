generate_signals

t_span = [-50,50];

sinc_filter = piecewise(t == 0, 1, sinc(t*0.7/pi))*0.7/pi;

%% Plot
fig = figure;
defaultPos = get(groot, 'DefaultFigurePosition'); % Get MATLAB's default figure size

% Modify height to be half of the default
newHeight = defaultPos(4) *0.7;

% Set new figure position with the same width but half the height
set(fig, 'Position', [defaultPos(1), defaultPos(2) + defaultPos(4) / 2, defaultPos(3), newHeight]);

hold on;
fplot(sinc_filter,t_span, 'r', 'LineWidth', 6, 'DisplayName', 'Sinc filter');
fplot(filter_signal,t_span, '-.b', 'LineWidth', 6, 'DisplayName', 'Cosine filter');

% Customize plot
xlabel(['Time [2' char(960) '/' char(969) ']'])
ylabel('Amplitude [arb. units]');
legend('show','fontweight','bold','fontsize',12);
grid on;
ax=gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
ylim([-0.052,0.2388])
hold off;
