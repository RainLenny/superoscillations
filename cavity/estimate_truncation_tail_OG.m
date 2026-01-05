function [eps_N, r_est, eps_bound] = estimate_truncation_tail_OG( ...
    omega, nu0, N, J_fluc, J_drive, Drive_integral, T_final)

% Estimate truncation tail assuming geometric decay
%
% Outputs:
%   eps_N     : ||Pe(N+1) - Pe(N)||_inf
%   r_est     : estimated geometric decay ratio
%   eps_bound : bound on ||Pe(inf) - Pe(N)||_inf

fprintf('\n--- Truncation tail estimation ---\n');

% Run N, N+1, N+2 (no recursion!)
[tN,  PeN,  ~] = multimode_JC_driven( ...
    omega, nu0, N,   J_fluc, J_drive, Drive_integral, T_final, false);

[tNp1,PeNp1,~] = multimode_JC_driven( ...
    omega, nu0, N+1, J_fluc, J_drive, Drive_integral, T_final, false);

[tNp2,PeNp2,~] = multimode_JC_driven( ...
    omega, nu0, N+2, J_fluc, J_drive, Drive_integral, T_final, false);

% Interpolate to common grid
PeNp1_i = interp1(tNp1, PeNp1, tN,   'linear');
PeNp2_i = interp1(tNp2, PeNp2, tNp1, 'linear');

% Truncation increments
eps_N   = max(abs(PeNp1_i - PeN));
eps_Np1 = max(abs(PeNp2_i - PeNp1));

% Decay ratio estimate
r_est = eps_Np1 / eps_N;

% Safety: avoid nonsense bounds
if r_est >= 1
    warning('Truncation increments not decaying (r >= 1). Bound invalid.');
    eps_bound = Inf;
else
    eps_bound = eps_N / (1 - r_est);
end

fprintf('eps_N      = %.3e\n', eps_N);
fprintf('eps_N+1    = %.3e\n', eps_Np1);
fprintf('r_est      = %.3f\n', r_est);
fprintf('BOUND inf  = %.3e\n', eps_bound);
fprintf('---------------------------------\n');

end
