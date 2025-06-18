superoscilation = ones(1,6);
w_axis_superoscilation = [0.4,0.45,0.5,0.55,0.6,0.65];

d_omega = 0.00001;  
w_axis = -2:d_omega:2;

sinc_signal = (w_axis >= -1) & (w_axis <= 1);  % Define sinc signal
sinc_signal = pi*sinc_signal;

%% Plot
fig = figure;
defaultPos = get(groot, 'DefaultFigurePosition'); % Get MATLAB's default figure size

% Modify height to be half of the default
newHeight = defaultPos(4) *0.6;

% Set new figure position with the same width but half the height
set(fig, 'Position', [defaultPos(1), defaultPos(2) + defaultPos(4) / 2, defaultPos(3), newHeight]);

hold on;

% Plot the dots for superoscilation
plot_superoscillations = plot(w_axis_superoscilation, superoscilation, 'o', 'color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

% Draw vertical lines from each dot to the x-axis
for i = 1:length(w_axis_superoscilation)
    line([w_axis_superoscilation(i), w_axis_superoscilation(i)], [0, superoscilation(i)], 'color', 'r', 'LineWidth', 2);
end

% Plot sinc signal
plot_sinc = plot(w_axis, sinc_signal, 'b-', 'LineWidth', 4);



% Customize plot
xlabel(['Frequency [' char(969) ']'])
ylabel('Amplitude [arb. units]');
legend([plot_superoscillations, plot_sinc], {'Superoscilation', 'Sinc'}, 'fontweight', 'bold', 'fontsize', 12);
grid on;
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
xlim([-1.5,1.5])
hold off;
