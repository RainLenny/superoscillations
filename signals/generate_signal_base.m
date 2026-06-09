function [signal, angular_freqs] = generate_signal_base(angular_freqs, amps, use_gaussian, use_conj)
%GENERATE_SIGNAL_BASE Base signal generation function using sum-of-exponentials logic
%
%   [signal, angular_freqs] = generate_signal_base(angular_freqs, amps, use_gaussian, use_conj)
%
%   INPUTS:
%     angular_freqs : row vector of angular frequencies
%     amps          : row vector of amplitudes (default 1)
%     use_gaussian  : boolean, whether to apply the Gaussian canvas (default true)
%     use_conj      : boolean, whether to keep complex conjugate signal (default false)
%
%   OUTPUTS:
%     signal        : function handle representing the generated signal
%     angular_freqs : 1xK vector of frequencies

    if nargin < 2 || isempty(amps)
        amps = 1;
    end
    if nargin < 3 || isempty(use_gaussian)
        use_gaussian = true;
    end
    if nargin < 4 || isempty(use_conj)
        use_conj = false;
    end

    angular_freqs = angular_freqs(:).';
    amps = amps(:).';

    t_0 = 250;
    T = 100;

    signal = @(t) reshape(conj(sum(amps .* exp(1i * angular_freqs .* t(:)), 2)), size(t));

    if use_gaussian
        signal = @(t) reshape(signal(t) .* exp(-(t(:)-t_0).^2./T^2), size(t));
    end
    if ~use_conj
        signal = @(t) real(signal(t));
    end

end
