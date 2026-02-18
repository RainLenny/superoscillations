%% CONSTANTS
clear; clc
freq_scaling = 1;
amp_scaling = 1;
[SO_signal, Cos_signal,angular_freqs_SO, angular_freqs_COS]  = generate_signals_Baranov_2014(freq_scaling,amp_scaling);

nu0 = 1;
n_max = 10;
factor = 1;
J_drive = 60/factor;
J_fluc = 0.001*factor;

T_final = 500;

%% Dyamics computation

% SO
[tgrid_SO, Pe_SO] = multimode_JC_driven(angular_freqs_SO, nu0, n_max, J_fluc, J_drive, SO_signal, T_final);
[tgrid_SO_no_cavity, Pe_SO_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, SO_signal, [0,T_final]);

%COS
[tgrid_COS, Pe_COS] = multimode_JC_driven(angular_freqs_SO, nu0, n_max, J_fluc, J_drive, Cos_signal, T_final);
[tgrid_COS_no_cavity, Pe_COS_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, Cos_signal, [0,T_final]);

%% Error computation


%% PLOT
figure;
ax1 = nexttile;
hold on;

plot(tgrid_SO, Pe_SO, 'LineWidth', 5, 'Color', 'red');
plot(tgrid_SO_no_cavity, Pe_SO_no_cavity, '--', 'LineWidth', 3, 'Color', 'k');

plot(tgrid_COS, Pe_COS, 'LineWidth', 5, 'Color', 'blue');
plot(tgrid_COS_no_cavity, Pe_COS_no_cavity, '-.', 'LineWidth', 3, 'Color', 'k');

grid on;
xlabel('Time [arb]');
ylabel('Excitation probability');
ax1.FontWeight = 'bold';
ax1.FontSize   = 12;
legend({'SO','SO no fluc','COS','COS no fluc'}, 'Location','best');

% plot_signals(freq_scaling,amp_scaling)