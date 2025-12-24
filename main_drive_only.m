%% CONSTANTS
clear; clc
freq_scaling = 1;
amp_scaling = 1;
[SO_signal, Cos_signal,angular_freqs_SO, angular_freqs_COS]  = generate_signals_Baranov_2014(freq_scaling,amp_scaling);

nu0 = 1;
tspan = [0,500];
Jtot = 0.03;


[tgrid_SO, Pe_SO] = JC_drive_only(nu0, Jtot, SO_signal, tspan);
[tgrid_COS, Pe_COS] = JC_drive_only(nu0, Jtot, Cos_signal, tspan);

%% PLOTS
figure;
ax1 = nexttile;
hold on;

plot(tgrid_SO, Pe_SO, 'LineWidth',5,'Color','red');
plot(tgrid_COS, Pe_COS, 'LineWidth',5,'Color','blue');

grid on;
xlabel('Time [arb]');
ylabel('Excited state probability');
ax1.FontWeight = 'bold';
ax1.FontSize   = 12;
legend({'SO','COS'}, 'Location','best');


plot_signals(freq_scaling,amp_scaling)