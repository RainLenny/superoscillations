function [Cos_signal, angular_freqs_COS] = generate_Cos_reference(angular_freqs_COS, amps_COS, varargin)
%GENERATE_COS_REFERENCE Generate a general-purpose Cosine reference signal
%
%   [Cos_signal, angular_freqs_COS] = generate_Cos_reference(angular_freqs_COS, amps_COS, varargin)
%
%   INPUTS:
%     angular_freqs_COS : scalar or vector of angular frequencies
%     amps_COS          : scalar or vector of amplitudes (default 1)
%     varargin          : optional arguments passed directly to generate_signal_base
%                         e.g., use_gaussian, use_conj, use_pulse_shaping
%
%   OUTPUTS:
%     Cos_signal        : function handle representing the Cosine reference signal
%     angular_freqs_COS : 1xK vector of frequencies

    if nargin < 2 || isempty(amps_COS)
        amps_COS = 1;
    end

    [Cos_signal, angular_freqs_COS] = generate_signal_base(angular_freqs_COS, amps_COS, varargin{:});

end
