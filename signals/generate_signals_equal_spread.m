function [SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS] = generate_signals_equal_spread(freq_scaling, amp_scaling, norm_method)

% Default to effective drive normalization if not explicitly provided
if nargin < 3
    norm_method = 'effective'; 
end

% Parameters
load('SO_signal_data.mat', 'amps_SO', 'angular_freqs_SO');

% Contain the same frequency components
angular_freqs_COS = angular_freqs_SO;

t_0 = 250;
T = 100;

% Function handles for fast-oscillating components
% Separated to easily compute the external drive's envelope derivative
s_t_SO = @(t) sum(amps_SO .* exp(1i * angular_freqs_SO .* t(:)), 2);
s_t_COS_unscaled = @(t) sum(1 .* exp(1i * angular_freqs_COS .* t(:)), 2);

% Function handles for Effective Drives (with Gaussian envelope)
SO_signal = @(t) reshape(conj(s_t_SO(t) .* exp(-(t(:)-t_0).^2./T^2)), size(t));
Cos_signal_unscaled = @(t) reshape(conj(s_t_COS_unscaled(t) .* exp(-(t(:)-t_0).^2./T^2)), size(t));

% --- NORMALIZATION LOGIC ---
% Create a dense time grid around the pulse to numerically find the signal peaks
t_grid = linspace(t_0 - 5*T, t_0 + 5*T, 100000)'; 

switch lower(norm_method)
    case 'effective'
        % Normalize by the peak of the Effective Drive f(t)
        max_eff_SO = max(abs(SO_signal(t_grid)));
        max_eff_COS_unscaled = max(abs(Cos_signal_unscaled(t_grid)));
        
        scale_factor = max_eff_SO / max_eff_COS_unscaled;
        target_peak_SO = max_eff_SO; % Stored for unit test
        
    case 'external'
        % Normalize by the peak of the External Drive D(t)
        % Envelope derivative: 2(t-t_0)/T^2 * exp(...)
        ext_env = @(t) (2 .* (t - t_0) ./ T^2) .* exp(-(t - t_0).^2 ./ T^2);
        
        % Calculate peaks of the physical external drives
        max_ext_SO = max(abs(ext_env(t_grid) .* s_t_SO(t_grid)));
        max_ext_COS_unscaled = max(abs(ext_env(t_grid) .* s_t_COS_unscaled(t_grid)));
        
        scale_factor = max_ext_SO / max_ext_COS_unscaled;
        target_peak_SO = max_ext_SO; % Stored for unit test
        
    otherwise
        error('Invalid norm_method. Please use ''effective'' or ''external''.');
end

% Apply scaling factor to the uniform amplitudes
amps_COS = 1 * scale_factor;

% Final scaled Cosine Signal (Effective Drive)
Cos_signal = @(t) reshape(conj(sum(amps_COS .* exp(1i * angular_freqs_COS .* t(:)), 2) .* exp(-(t(:)-t_0).^2./T^2)), size(t));


%% INBUILT VERIFICATION (UNIT TEST)
% Check if the peaks match within a tight numerical tolerance depending on the chosen method

tolerance = 1e-8;

if strcmpi(norm_method, 'effective')
    peak_Cos_final = max(abs(Cos_signal(t_grid)));
else
    % Calculate the final external drive peak for the scaled Cos signal
    s_t_COS_final = @(t) sum(amps_COS .* exp(1i * angular_freqs_COS .* t(:)), 2);
    peak_Cos_final = max(abs(ext_env(t_grid) .* s_t_COS_final(t_grid)));
end

if abs(target_peak_SO - peak_Cos_final) > tolerance
    % Throw a custom error breaking execution if the math failed
    error('SignalGen:PeakMismatch', ...
        'Validation Failed: Peaks do not match for %s normalization.\nPeak SO: %.10f\nPeak COS: %.10f\nDifference: %.10f', ...
        norm_method, target_peak_SO, peak_Cos_final, abs(target_peak_SO - peak_Cos_final));
end

end