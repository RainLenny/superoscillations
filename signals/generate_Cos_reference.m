function [Cos_signal, angular_freqs_COS, amps_COS] = generate_Cos_reference(angular_freqs_COS, amps_COS, use_gaussian, use_conj)
%GENERATE_COS_REFERENCE Generate a general-purpose Cosine reference signal
%
%   [Cos_signal, angular_freqs_COS, amps_COS] = generate_Cos_reference(angular_freqs_COS, amps_COS, use_gaussian, use_conj)
%
%   INPUTS:
%     angular_freqs_COS : scalar or vector of angular frequencies
%     amps_COS          : scalar or vector of amplitudes (default 1)
%     use_gaussian      : boolean, whether to apply the Gaussian canvas (default true)
%     use_conj          : boolean, whether to keep complex conjugate signal (default false)
%
%   OUTPUTS:
%     Cos_signal        : function handle representing the Cosine reference signal
%     angular_freqs_COS : 1xK vector of frequencies
%     amps_COS          : 1xK vector of complex amplitudes

    if nargin < 2 || isempty(amps_COS)
        amps_COS = 1;
    end
    if nargin < 3 || isempty(use_gaussian)
        use_gaussian = true;
    end
    if nargin < 4 || isempty(use_conj)
        use_conj = false;
    end

    [Cos_signal, angular_freqs_COS, amps_COS] = generate_signal_base(angular_freqs_COS, amps_COS, use_gaussian, use_conj);

end
