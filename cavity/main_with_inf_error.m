%% CONSTANTS
clear; clc
freq_scaling = 1;
amp_scaling = 1;
[SO_signal, Cos_signal,angular_freqs_SO, angular_freqs_COS]  = generate_signals_Baranov_2014(freq_scaling,amp_scaling);

nu0 = 1;
nmax = 3;
J_drive = 30;
J_fluc = 0.001;

T_final = 600;
do_err_est = 1;

%% Dyamics computation

% SO
[tgrid_SO, Pe_SO,eps_trunc_SO] = multimode_JC_driven(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, SO_signal, T_final, do_err_est);
[tgrid_SO_no_cavity, Pe_SO_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, SO_signal, [0,T_final]);

%COS
[tgrid_COS, Pe_COS,eps_trunc_COS] = multimode_JC_driven(angular_freqs_SO, nu0, nmax, J_fluc, J_drive, Cos_signal, T_final, do_err_est);
[tgrid_COS_no_cavity, Pe_COS_no_cavity] = JC_drive_only(nu0, J_fluc*J_drive, Cos_signal, [0,T_final]);


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

%% INFINITE-HILBERT-SPACE ERROR BOUND

N = nmax;   % current truncation
[epsN_SO, r_SO, epsInf_SO] = estimate_truncation_tail_OG( ...
    angular_freqs_SO, nu0, N, J_fluc, J_drive, SO_signal, T_final);

[epsN_COS, r_COS, epsInf_COS] = estimate_truncation_tail_OG( ...
    angular_freqs_COS, nu0, N, J_fluc, J_drive, Cos_signal, T_final);

fprintf('\n===== INFINITE-DIMENSION ERROR BOUNDS =====\n');
fprintf('SO   : ||Pe(N)-Pe(inf)|| <= %.3e   (r = %.3f)\n', epsInf_SO,  r_SO);
fprintf('COS  : ||Pe(N)-Pe(inf)|| <= %.3e   (r = %.3f)\n', epsInf_COS, r_COS);
fprintf('==========================================\n');
