function normalized_signals = normalize_signal_list_by_energy(signals_list, t_scan)
    % NORMALIZE_SIGNAL_LIST Normalizes a list of signal function handles.
    % This version uses arrayfun to safely handle function handles that 
    % might not support vector inputs (avoiding dimension mismatch errors).

    normalized_signals = signals_list; 
    
    if isempty(signals_list)
        return;
    end

    % --- Helper to safely evaluate function handles ---
    % This forces element-by-element evaluation, preventing 
    % "Incompatible sizes" errors if the handle isn't vectorized.
    evaluate_safe = @(func, t) arrayfun(func, t);

    % 1. Calculate the Energy of the Reference Signal (Signal 1)
    ref_func = signals_list{1};
    
    % Use the safe evaluator
    ref_vals = evaluate_safe(ref_func, t_scan);
    
    % Calculate Energy (L2 Norm)
    ref_energy = sqrt(sum(abs(ref_vals).^2));
    
    fprintf('Reference Signal (Signal 1) Energy Metric: %.4f\n', ref_energy);

    % 2. Loop through the remaining signals and normalize
    for i = 2:length(signals_list)
        curr_func = signals_list{i};
        
        % Use the safe evaluator
        curr_vals = evaluate_safe(curr_func, t_scan);
        
        curr_energy = sqrt(sum(abs(curr_vals).^2));
        
        if curr_energy == 0
            scale_factor = 1; 
        else
            scale_factor = ref_energy / curr_energy;
        end
        
        % Update the handle with the scaling factor
        normalized_signals{i} = @(t) curr_func(t) * scale_factor;
        
        fprintf('Signal %d normalized by factor: %.4f\n', i, scale_factor);
    end
end