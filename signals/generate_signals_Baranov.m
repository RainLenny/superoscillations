function [SO_signal, Cos_signal, angular_freqs_SO, angular_freqs_COS] = generate_signals_Baranov(freq_scaling, amp_scaling, norm_method)

% Default to effective drive normalization if not explicitly provided
if nargin < 3
    norm_method = 'effective'; 
end

% Parameters
load('SO_signal_data.mat', 'amps_SO', 'angular_freqs_SO');

angular_freqs_COS = 0.7 * freq_scaling;

t_0 = 250;
T = 100;

% Function handle for superoscillatory signal (Effective Drive)
% (Note: works best when t is a column vector due to sum(..., 2))
SO_signal = @(t) conj(sum( amps_SO .* exp(1i * angular_freqs_SO .* t), 2) .* exp(-(t-t_0).^2./T^2));

% --- NORMALIZATION LOGIC ---
% Create a dense time grid around the pulse to numerically find the signal peaks
t_grid = linspace(t_0 - 5*T, t_0 + 5*T, 100000)'; 

switch lower(norm_method)
    case 'effective'
        % 1. Normalize by the Effective Drive f(t)
        % We equate the peak amplitude of the SO effective drive to the COS effective drive.
        
        SO_eff_eval = SO_signal(t_grid);
        max_eff_SO = max(abs(SO_eff_eval));
        
        % The Gaussian envelope of Cos_signal peaks at 1, so its max amplitude is exactly Cos_normalization
        Cos_normalization = max_eff_SO;
        
    case 'external'
        % 2. Normalize by the External Drive D(t)
        % Based on Eq. 8 in your paper, the external drive is the derivative of the Gaussian 
        % envelope multiplied by the fast-oscillating components.
        
        % External envelope factor: 2(t-t_0)/T^2 * exp(...)
        ext_env = @(t) (2 .* (t - t_0) ./ T^2) .* exp(-(t - t_0).^2 ./ T^2);
        
        % The fast-oscillating part of the SO signal
        s_t_SO = @(t) sum(amps_SO .* exp(1i * angular_freqs_SO .* t), 2);
        
        % Calculate the peak of the physical external SO drive
        D_SO_mag = abs(ext_env(t_grid) .* s_t_SO(t_grid));
        max_ext_SO = max(D_SO_mag);
        
        % Calculate the peak of the physical external COS drive's envelope
        max_ext_env = max(abs(ext_env(t_grid)));
        
        % Equate the peaks: Cos_normalization * max_ext_env = max_ext_SO
        Cos_normalization = max_ext_SO / max_ext_env;
        
    otherwise
        error('Invalid norm_method. Please use ''effective'' or ''external''.');
end

% Function handle for the reference Cosine signal (Effective Drive)
Cos_signal = @(t) -conj(Cos_normalization .* exp(1i * angular_freqs_COS .* t) .* exp(-(t-t_0).^2./T^2));

end