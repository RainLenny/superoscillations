function [varargout] = normalize_signals(signals, norm_type, t_guess)
%NORMALIZE_SIGNALS Normalize a list of signals (Real or Complex) to match the first signal.
%
%   [S1_norm, S2_norm, ...] = normalize_signals(signals)
%   [S1_norm, S2_norm, ..., Rel_N] = normalize_signals(signals, norm_type)
%   [..., Rel_N] = normalize_signals(signals, norm_type, t_guess)
%   normalized_cell = normalize_signals(signals, ...)
%   [S_norm, Rel_N] = normalize_signals(S, ...)
%
%   Supported norm_type values:
%       'peak'   - Normalizes based on peak amplitude (default)
%       'energy' - Normalizes based on signal energy
%       'pass'   - Does not change the signals (returns original signals, Rel_N = 1)

%% 1. Input Parsing & Standardization
if nargin < 2 || isempty(norm_type), norm_type = 'peak'; end
if nargin < 3 || isempty(t_guess),   t_guess = [-1000, 1000]; end

% Set search bounds
if isscalar(t_guess)
    t_search = [t_guess - 1000, t_guess + 1000];
else
    t_search = [t_guess(1), t_guess(end)];
end

% Standardize input to cell array
was_single_handle = isa(signals, 'function_handle');
if was_single_handle
    signals = {signals};
elseif ~iscell(signals)
    error('Input ''signals'' must be a cell array of function handles or a single function handle.');
end

num_signals = length(signals);
Raw_N = ones(1, num_signals); % Default to 1 to gracefully handle failures and 'pass' mode

%% 2. Pre-allocation & Strategy Setup
norm_type = lower(norm_type);

if strcmp(norm_type, 'peak')
    options = optimset('Display', 'off', 'TolX', 1e-10);
    % Pre-generate grid once outside the loop to save overhead
    t_samples = linspace(t_search(1), t_search(2), 100000);
elseif ~strcmp(norm_type, 'energy') && ~strcmp(norm_type, 'pass')
    error('Unknown norm_type: ''%s''. Please use ''peak'', ''energy'', or ''pass''.', norm_type);
end

%% 3. Normalization Calculation
for i = 1:num_signals
    S = signals{i};

    if strcmp(norm_type, 'pass')
        % Keep Raw_N(i) = 1 to pass the signal through unchanged
        continue;
    elseif strcmp(norm_type, 'peak')
        % 1. Global search via dense sampling
        try
            samp_vals = abs(S(t_samples));
        catch
            samp_vals = arrayfun(@(t) abs(S(t)), t_samples);
        end

        [max_val, max_idx] = max(samp_vals);

        if max_val > 0
            % Define a tight bounding box around the sampled peak
            t_lower = t_samples(max(1, max_idx - 2));
            t_upper = t_samples(min(end, max_idx + 2));

            % 2. Exact mathematical refinement using bounded 1D optimization
            [~, min_val] = fminbnd(@(t) -abs(S(t)), t_lower, t_upper, options);
            Raw_N(i) = abs(min_val);
        else
            warning('Signal %d has a peak of 0 or invalid values. Raw factor set to 1.', i);
        end

    else % 'energy'
        try
            E = integral(@(t) abs(S(t)).^2, t_search(1), t_search(2), ...
                'ArrayValued', true, 'RelTol', 1e-8, 'AbsTol', 1e-10);

            if E > 0
                Raw_N(i) = sqrt(E);
            else
                warning('Signal %d has zero or negative energy. Raw factor set to 1.', i);
            end
        catch ME
            warning('Integration failed for signal %d: %s. Raw factor set to 1.', i, ME.message);
        end
    end
end

%% 4. Relative Scaling & Handle Generation
% Compute relative factors against the reference signal
ref_factor = Raw_N(1);
Rel_N = Raw_N ./ ref_factor;

S_norm_cell = cell(1, num_signals);
for i = 1:num_signals
    % Securely isolate both the specific signal and its scale factor
    S_curr = signals{i};
    Rel_curr = Rel_N(i);
    S_norm_cell{i} = @(t) S_curr(t) ./ Rel_curr;
end

%% 5. Output Mapping (varargout)
if was_single_handle
    varargout{1} = S_norm_cell{1};
    if nargout > 1, varargout{2} = Rel_N; end
else
    if nargout == num_signals || nargout == num_signals + 1
        varargout(1:num_signals) = S_norm_cell;
        if nargout == num_signals + 1
            varargout{num_signals + 1} = Rel_N;
        end
    else
        varargout{1} = S_norm_cell;
        if nargout > 1, varargout{2} = Rel_N; end
    end
end
end