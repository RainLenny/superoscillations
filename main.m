%% CONSTANTS
clear; clc
freq_scaling = 1;
amp_scaling = 1;
[SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS]  = generate_signals_Baranov_2014(freq_scaling,amp_scaling);

nu0 = 1;
nmax = 0;
J_drive = 30;
J_fluc = 0.001;
T1 = 500;
T2 = 300;

T_final = 600;
do_err_est = 0;

%% Dyamics computation

% SO
[tgrid_SO, Pe_SO,eps_trunc_SO] = multimode_JC_driven(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, SO_signal, T_final, do_err_est);
% [tgrid_SO, Pe_SO,eps_trunc_SO] = master_multimode_JC_driven(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, SO_signal, T_final, T1, T2, do_err_est);
[tgrid_SO_no_cavity, Pe_SO_no_cavity] = master_JC_drive_only(nu0, J_fluc*J_drive, SO_signal, [0,T_final], T1, T2,0);

%COS
[tgrid_COS, Pe_COS,eps_trunc_COS] = multimode_JC_driven(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, Cos_signal, T_final, do_err_est);
% [tgrid_COS, Pe_COS,eps_trunc_COS] = master_multimode_JC_driven(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, Cos_signal, T_final, T1, T2, do_err_est);
[tgrid_COS_no_cavity, Pe_COS_no_cavity] = master_JC_drive_only(nu0, J_fluc*J_drive, Cos_signal, [0,T_final],T1, T2,0);

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
ylabel('Excited state probability');
ax1.FontWeight = 'bold';
ax1.FontSize   = 12;
legend({'SO','SO no fluc','COS','COS no fluc'}, 'Location','best');

fprintf('\n===== Numerical Errors =====\n');
fprintf('SO error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_SO);
fprintf('COS error (Hilbert truncation nmax+1): %.8e\n', eps_trunc_COS);
fprintf('============================\n');

% plot_signals(freq_scaling,amp_scaling)