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
%     t_guess   : (Optional) Initial guess for peak search if exact 
%                 analytical solving fails (default: 250).
%
%   OUTPUTS:
%     normalized function handles, plus the relative normalization factor Rel_N.
%     (Dividing a signal by its Rel_N makes it match the first signal).

    %% 1. Input Parsing
    if nargin < 2 || isempty(norm_type)
        norm_type = 'peak';
    end
    if nargin < 3 || isempty(t_guess)
        t_guess = 250;
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
    Raw_N = zeros(1, num_signals); % Stores the absolute peak/energy
    Rel_N = zeros(1, num_signals); % Stores the relative factor compared to signal 1
    S_norm_cell = cell(1, num_signals);

    %% 2. Normalization Calculation
    for i = 1:num_signals
        S = signals{i};
        
        switch lower(norm_type)
            case 'peak'
                % Use fminsearch to find the minimum of the negative absolute value
                obj_fun = @(t) -abs(S(t));
                options = optimset('Display', 'off'); % Suppress console output
                
                [~, min_val] = fminsearch(obj_fun, t_guess, options);
                peak_val = abs(min_val);
                
                if peak_val == 0
                    warning('Signal %d has a peak of 0. Raw factor set to 1 to prevent division by zero.', i);
                    Raw_N(i) = 1;
                else
                    Raw_N(i) = peak_val;
                end
                
            case 'energy'
                % Integrate |S(t)|^2 from -Inf to Inf
                % ArrayValued=true ensures it works even if S(t) isn't strictly vectorized
                energy_fun = @(t) abs(S(t)).^2;
                
                try
                    E = integral(energy_fun, -Inf, Inf, 'ArrayValued', true);
                catch ME
                    warning('Integration failed for signal %d: %s. Raw factor set to 1.', i, ME.message);
                    E = 0;
                end
                
                if E == 0
                    warning('Signal %d has zero energy. Raw factor set to 1.', i);
                    Raw_N(i) = 1;
                else
                    Raw_N(i) = sqrt(E); % Amplitude norm factor is the square root of energy
                end
                
            otherwise
                error('Unknown norm_type: ''%s''. Please use ''peak'' or ''energy''.', norm_type);
        end
        
        % Set the reference factor on the first iteration
        if i == 1
            ref_factor = Raw_N(1);
        end
        
        % Calculate relative scaling factor (Targeting the first signal's amplitude/energy)
        Rel_N(i) = Raw_N(i) / ref_factor;
        
        % Create the normalized function handle.
        % We capture Rel_N(i) in a local scalar variable (Rel_curr) so the 
        % anonymous function closes over the scalar correctly.
        Rel_curr = Rel_N(i);
        S_norm_cell{i} = @(t) S(t) ./ Rel_curr;
    end

    %% 3. Output Formatting (varargout)
    if was_single_handle
        % Return single handle + factor
        varargout{1} = S_norm_cell{1};
        if nargout > 1
            varargout{2} = Rel_N;
        end
    else
        % If nargout matches the exact number of signals (plus optional Rel_N)
        if nargout == num_signals || nargout == num_signals + 1
            for i = 1:num_signals
                varargout{i} = S_norm_cell{i};
            end
            if nargout == num_signals + 1
                varargout{num_signals + 1} = Rel_N;
            end
        else
            % Default fallback: return the cell array + factor
            varargout{1} = S_norm_cell;
            if nargout > 1
                varargout{2} = Rel_N;
            end
        end
    end
end