function [SO_signal, angular_freqs_SO] = generate_SO_leakage_comparison(freq_scaling, amp_scaling, varargin)
%GENERATE_SO_LEAKAGE_COMPARISON Generate the Superoscillating (SO) signal for Leakage Comparison
%
%   [SO_signal, angular_freqs_SO] = generate_SO_leakage_comparison(freq_scaling, amp_scaling, varargin)
%
%   INPUTS:
%     freq_scaling  : multiplier for frequencies (default 1)
%     amp_scaling   : multiplier for amplitudes (default 1)
%     varargin      : optional arguments passed directly to generate_signal_base
%                     e.g., use_gaussian, use_conj, use_pulse_shaping
%
%   OUTPUTS:
%     SO_signal         : function handle of the SO signal
%     angular_freqs_SO  : 1xN vector of scaled angular frequencies

    if nargin < 1 || isempty(freq_scaling)
        freq_scaling = 1;
    end
    if nargin < 2 || isempty(amp_scaling)
        amp_scaling = 1;
    end
    
    % Maintain original default behavior: use_gaussian=true, use_conj=true
    if length(varargin) < 1 || isempty(varargin{1})
        varargin{1} = true; % use_gaussian
    end
    if length(varargin) < 2 || isempty(varargin{2})
        varargin{2} = true; % use_conj
    end

    % Import from the .mat file instead of hardcoding
    data = load('SO_Baranov.mat', 'amps_SO', 'angular_freqs_SO');
    amps_SO = data.amps_SO * amp_scaling;
    angular_freqs_SO = data.angular_freqs_SO * freq_scaling;

    [SO_signal, angular_freqs_SO] = generate_signal_base(angular_freqs_SO, amps_SO, varargin{:});

end
