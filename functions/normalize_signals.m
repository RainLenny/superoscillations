function [varargout] = normalize_signals(signals, norm_type, t_guess)
%NORMALIZE_SIGNALS Normalize a list of signals (Real or Complex) to match the first signal.
%
%   [S1_norm, S2_norm, ...] = normalize_signals(signals)
%   [S1_norm, S2_norm, ..., Rel_N] = normalize_signals(signals, norm_type)
%   [..., Rel_N] = normalize_signals(signals, norm_type, t_guess)
%   normalized_cell = normalize_signals(signals, ...)
%   [S_norm, Rel_N] = normalize_signals(S, ...)
%
%   INPUTS:
%     signals   : cell array of function handles OR a single function handle.
%                 The FIRST signal in this list is the reference signal.
%     norm_type : 'peak' (default) or 'energy'
%     t_guess   : (Optional) Search space. Can be a scalar guess 
%                 (searches +/- 1000 around it) or a 2-element vector [t_min, t_max].
%                 (default: [-1000, 1000]).
%
%   OUTPUTS:
%     normalized function handles, plus the relative normalization factor Rel_N.

    %% 1. Input Parsing
    if nargin < 2 || isempty(norm_type)
        norm_type = 'peak';
    end
    
    % Determine the search bounds for peak finding/integration
    if nargin < 3 || isempty(t_guess)
        t_search = [-1000, 1000];
    elseif isscalar(t_guess)
        t_search = [t_guess - 1000, t_guess + 1000];
    elseif length(t_guess) >= 2
        t_search = [t_guess(1), t_guess(end)];
    end

    % Standardize input to a cell array
    was_single_handle = false;
    if isa(signals, 'function_handle')
        was_single_handle = true;
        signals = {signals};
    elseif ~iscell(signals)
        error('Input ''signals'' must be a cell array of function handles or a single function handle.');
    end

    num_signals = length(signals);
    Raw_N = zeros(1, num_signals); 
    Rel_N = zeros(1, num_signals); 
    S_norm_cell = cell(1, num_signals);

    %% 2. Pre-allocate Options for Speed
    if strcmpi(norm_type, 'peak')
        % fminbnd is much faster/more accurate for 1D than fminsearch
        options = optimset('Display', 'off', 'TolX', 1e-10);
        % Higher density grid to prevent stepping over narrow oscillations
        N_samples = 100000;
        t_samples = linspace(t_search(1), t_search(2), N_samples);
    end

    %% 3. Normalization Calculation
    for i = 1:num_signals
        S = signals{i};
        
        switch lower(norm_type)
            case 'peak'
                % 1. Global search via dense sampling
                try
                    samp_vals = abs(S(t_samples));
                catch
                    samp_vals = arrayfun(@(t) abs(S(t)), t_samples);
                end
                
                [~, max_idx] = max(samp_vals);
                
                % Define a tight bounding box around the sampled peak
                idx_min = max(1, max_idx - 2);
                idx_max = min(N_samples, max_idx + 2);
                t_lower = t_samples(idx_min);
                t_upper = t_samples(idx_max);
                
                % 2. Exact mathematical refinement using bounded 1D optimization
                obj_fun = @(t) -abs(S(t));
                [~, min_val] = fminbnd(obj_fun, t_lower, t_upper, options);
                peak_val = abs(min_val);
                
                if peak_val == 0
                    warning('Signal %d has a peak of 0. Raw factor set to 1.', i);
                    Raw_N(i) = 1;
                else
                    Raw_N(i) = peak_val;
                end
                
            case 'energy'
                energy_fun = @(t) abs(S(t)).^2;
                
                try
                    % Bounded integration with tight tolerances to handle rapid oscillations
                    E = integral(energy_fun, t_search(1), t_search(2), ...
                        'ArrayValued', true, 'RelTol', 1e-8, 'AbsTol', 1e-10);
                catch ME
                    warning('Integration failed for signal %d: %s. Raw factor set to 1.', i, ME.message);
                    E = 0;
                end
                
                if E == 0
                    warning('Signal %d has zero energy. Raw factor set to 1.', i);
                    Raw_N(i) = 1;
                else
                    Raw_N(i) = sqrt(E);
                end
                
            otherwise
                error('Unknown norm_type: ''%s''. Please use ''peak'' or ''energy''.', norm_type);
        end
        
        % Set the reference factor on the first iteration
        if i == 1
            ref_factor = Raw_N(1);
        end
        
        % Calculate relative scaling factor
        Rel_N(i) = Raw_N(i) / ref_factor;
        
        % Capture Rel_N(i) locally to avoid closure issues in anonymous functions
        Rel_curr = Rel_N(i);
        S_norm_cell{i} = @(t) S(t) ./ Rel_curr;
    end

    %% 4. Output Formatting (varargout)
    if was_single_handle
        varargout{1} = S_norm_cell{1};
        if nargout > 1
            varargout{2} = Rel_N;
        end
    else
        if nargout == num_signals || nargout == num_signals + 1
            for i = 1:num_signals
                varargout{i} = S_norm_cell{i};
            end
            if nargout == num_signals + 1
                varargout{num_signals + 1} = Rel_N;
            end
        else
            varargout{1} = S_norm_cell;
            if nargout > 1
                varargout{2} = Rel_N;
            end
        end
    end
end