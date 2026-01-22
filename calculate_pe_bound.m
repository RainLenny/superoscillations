function bound = calculate_pe_bound(K, J, t, n_max)
% CALCULATE_PE_BOUND Computes the upper bound for Delta Pe(t)
%
%   Formula: 2 * sum_{n=n_max+1}^{inf} [ (sqrt(K)*J*t)^n / sqrt(n!) ]
%   Method: Uses Log-Sum-Exp trick to ensure numerical stability for large n.

    X = sqrt(K) * J * t;
    
    running_sum = 0;
    n = n_max + 1;
    
    % Stop when additions are smaller than machine precision relative to sum
    tolerance = 1e-16; 
    
    while true
        % Compute log(term) first to avoid factorial overflow
        % term = X^n / sqrt(n!)
        log_term = n * log(X) - 0.5 * gammaln(n + 1);
        
        term = exp(log_term);
        
        % Accumulate
        old_sum = running_sum;
        running_sum = running_sum + term;
        
        % Convergence check: If adding the term didn't change the sum significantly
        if (running_sum > 0 && (running_sum - old_sum) / running_sum < tolerance)
            break;
        elseif running_sum == 0 && term < tolerance
            break;
        end
        
        % Safety break
        if n > n_max + 2000
            break; 
        end
        
        n = n + 1;
    end
    
    bound = 2 * running_sum;
end