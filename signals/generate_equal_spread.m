function [SO_signal, angular_freqs_SO] = generate_equal_spread(input_freqs, freq_scaling, amp_scaling, varargin)
%GENERATE_SO_EQUAL_SPREAD Generate the Superoscillating (SO) signal for Equal Spread method
%
%   [SO_signal, angular_freqs_SO] = generate_SO_equal_spread(freq_scaling, amp_scaling, varargin)
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

    if nargin < 2 || isempty(freq_scaling)
        freq_scaling = 1;
    end
    if nargin < 3 || isempty(amp_scaling)
        amp_scaling = 1;
    end

    % Import from the .mat file instead of hardcoding
    amps_SO = amp_scaling;
    angular_freqs_SO = input_freqs * freq_scaling;

    [SO_signal, angular_freqs_SO] = generate_signal_base(angular_freqs_SO, amps_SO, varargin{:});

end
