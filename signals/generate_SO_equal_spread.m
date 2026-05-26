function [SO_signal, angular_freqs_SO] = generate_SO_equal_spread(freq_scaling, amp_scaling, use_gaussian)
%GENERATE_SO_EQUAL_SPREAD Generate the Superoscillating (SO) signal for Equal Spread method
%
%   [SO_signal, angular_freqs_SO] = generate_SO_equal_spread(freq_scaling, amp_scaling, use_gaussian)
%
%   INPUTS:
%     freq_scaling  : multiplier for frequencies (default 1)
%     amp_scaling   : multiplier for amplitudes (default 1)
%     use_gaussian  : boolean, whether to apply the Gaussian canvas (default true)
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

    % Import from the .mat file instead of hardcoding
    data = load('SO_Baranov.mat', 'amps_SO', 'angular_freqs_SO');
    amps_SO = amp_scaling;
    angular_freqs_SO = data.angular_freqs_SO * freq_scaling;

    t_0 = 250;
    T = 100;

    s_t_SO = @(t) sum(amps_SO .* exp(1i * angular_freqs_SO .* t(:)), 2);

    % Create regular or Gaussian canvas version
    if use_gaussian
        SO_signal = @(t) reshape(conj(s_t_SO(t) .* exp(-(t(:)-t_0).^2./T^2)), size(t));
    else
        SO_signal = @(t) reshape(conj(s_t_SO(t)), size(t));
    end
end
