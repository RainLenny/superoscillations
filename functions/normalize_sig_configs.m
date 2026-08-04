function sig_configs = normalize_sig_configs(sig_configs, method)
    % NORMALIZE_SIG_CONFIGS Normalizes the 'data' field in an array of sig_configs structs
    %
    % Usage:
    %   sig_configs = normalize_sig_configs(sig_configs, method)
    %
    % Inputs:
    %   sig_configs - struct array with field 'data'
    %   method      - (optional) 'energy' (default) or 'peak'
    %
    % Outputs:
    %   sig_configs - updated struct array with normalized 'data'
    
    if nargin < 2
        method = 'energy';
    end
    
    sigs = {sig_configs.data};
    [norm_sigs{1:length(sigs)}] = normalize_signals(sigs, method);
    
    for i = 1:length(sig_configs)
        sig_configs(i).data = norm_sigs{i};
    end
end
