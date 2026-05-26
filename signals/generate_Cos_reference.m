function [Cos_signal, angular_freqs_COS] = generate_Cos_reference(angular_freqs_COS, amps_COS, use_gaussian)
%GENERATE_COS_REFERENCE Generate a general-purpose Cosine reference signal
%
%   [Cos_signal, angular_freqs_COS] = generate_Cos_reference(angular_freqs_COS, amps_COS, use_gaussian)
%
%   INPUTS:
%     angular_freqs_COS : scalar or vector of angular frequencies
%     amps_COS          : scalar or vector of amplitudes (default 1)
%     use_gaussian      : boolean, whether to apply the Gaussian canvas (default true)
%
%   OUTPUTS:
%     Cos_signal        : function handle representing the Cosine reference signal
%     angular_freqs_COS : 1xK vector of frequencies

    if nargin < 2 || isempty(amps_COS)
        amps_COS = 1;
    end
    if nargin < 3 || isempty(use_gaussian)
        use_gaussian = true;
    end

    angular_freqs_COS = angular_freqs_COS(:).';
    amps_COS = amps_COS(:).';

    t_0 = 250;
    T = 100;

    % Create regular or Gaussian canvas version
    if use_gaussian
        Cos_signal = @(t) reshape(conj(sum(amps_COS .* exp(1i * angular_freqs_COS .* t(:)), 2) .* exp(-(t(:)-t_0).^2./T^2)), size(t));
    else
        Cos_signal = @(t) reshape(conj(sum(amps_COS .* exp(1i * angular_freqs_COS .* t(:)), 2)), size(t));
    end
end
