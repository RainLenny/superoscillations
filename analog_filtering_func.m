function analog_filtering_func(sym_signal,sym_superoscillation,angular_freqs, sym_filter, f_sampling,window_length_for_filter)
% All the syms signals given to the function should be in their time
% representation

%% Filtering
syms t real
filtered_signal = sym_conv(sym_signal,sym_filter,t);
filtered_superoscillation = sym_conv(sym_superoscillation,sym_filter,t);

%% Plots
%% Parameters
t_min = -15;    % start of time-axis
t_max =  15;    % end of time-axis
syms omega real
omega_lim =  2;     % ± limit for frequency axis (rad/s)

%% Plot time-domain signals
figure;
hold on;
fplot( real(sym_superoscillation), [t_min, t_max], '-','LineWidth',4,'Color', 'b', 'DisplayName','Superoscillation');
fplot(sym_signal, [t_min, t_max], '-','LineWidth',4,'Color', 'r', 'DisplayName','Signal');
fplot(abs(sym_filter), [t_min, t_max], '-','LineWidth',4,'Color', 'k', 'DisplayName','Filter');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Location','best');
grid on;

%% Plot symbolic FFTs (Fourier transforms)
% compute symbolic Fourier transforms
F_super = fourier(sym_superoscillation, t, omega);
F_signal  = fourier(sym_signal,    t, omega);
F_filt  = fourier(sym_filter,     t, omega);

figure;
hold on;
fplot(abs(F_signal),  [-omega_lim, omega_lim], '-','LineWidth',2,'Color', 'r', 'DisplayName','Signal');
fplot(abs(F_super),[-omega_lim, omega_lim], '-','LineWidth',2,'Color', 'b', 'DisplayName','Superoscillation');
fplot(abs(F_filt), [-omega_lim, omega_lim], '-','LineWidth',2,'Color', 'k', 'DisplayName','Filter');
xlabel('Frequency [\omega]');
ylabel('Magnitude');
title('Symbolic Fourier Transforms');
legend('Location','best');
grid on;
xlim([-omega_lim, omega_lim]);

%% Plot original vs filtered in tiled layout
figure;
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding     = 'compact';

% Superoscillating
ax1 = nexttile;
hold on;
fplot(real(sym_superoscillation),           [t_min, t_max], 'LineWidth',6,'Color', 'r', 'DisplayName','Original');
fplot(real(filtered_superoscillation),  [t_min, t_max], '--','LineWidth',6,'Color', 'b', 'DisplayName','Filtered');
title('Superoscillations');
legend('Location','best','FontWeight','bold','FontSize',12);
grid on;
xlim([t_min, t_max]);
ylim([-0.5,1]);
ax1.FontWeight = 'bold';
ax1.FontSize   = 12;

% Sinc
ax2 = nexttile;
hold on;
fplot(sym_signal,        [t_min, t_max], 'LineWidth',6,'Color', 'r', 'DisplayName','Original');
fplot(real(filtered_signal),      [t_min, t_max], '--','LineWidth',6,'Color', 'b', 'DisplayName','Filtered');
title('Signal');
legend('Location','best','FontWeight','bold','FontSize',12);
grid on;
xlim([t_min, t_max]);
ylim([-0.5,1]);
ax2.FontWeight = 'bold';
ax2.FontSize   = 12;

% Shared labels
xlabel(t, ['Time [2' char(960) '/' char(969) ']'],   'FontWeight','bold','FontSize',12);
ylabel(t, 'Amplitude [arb. units]',                  'FontWeight','bold','FontSize',12);

end

