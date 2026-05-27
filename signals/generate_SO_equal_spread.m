function [SO_signal, angular_freqs_SO] = generate_SO_equal_spread(freq_scaling, amp_scaling, use_gaussian, use_conj)
%GENERATE_SO_EQUAL_SPREAD Generate the Superoscillating (SO) signal for Equal Spread method
%
%   [SO_signal, angular_freqs_SO] = generate_SO_equal_spread(freq_scaling, amp_scaling, use_gaussian, use_conj)
%
%   INPUTS:
%     freq_scaling  : multiplier for frequencies (default 1)
%     amp_scaling   : multiplier for amplitudes (default 1)
%     use_gaussian  : boolean, whether to apply the Gaussian canvas (default true)
%     use_conj      : boolean, whether to keep complex conjugate signal (default true)
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
    if nargin < 3 || isempty(use_gaussian)
        use_gaussian = true;
    end
    if nargin < 4 || isempty(use_conj)
        use_conj = false;
    end

    % Import from the .mat file instead of hardcoding
    data = load('SO_Baranov.mat', 'amps_SO', 'angular_freqs_SO');
    amps_SO = amp_scaling;
    angular_freqs_SO = data.angular_freqs_SO * freq_scaling;

    [SO_signal, angular_freqs_SO] = generate_signal_base(angular_freqs_SO, amps_SO, use_gaussian, use_conj);

end
