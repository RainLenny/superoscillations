function [signal, angular_freqs, amps] = generate_signal_base(angular_freqs, amps, use_gaussian, use_conj, use_pulse_shaping)
%GENERATE_SIGNAL_BASE Base signal generation function using sum-of-exponentials logic
%
%   [signal, angular_freqs] = generate_signal_base(angular_freqs, amps, use_gaussian, use_conj, use_pulse_shaping)
%
%   INPUTS:
%     angular_freqs     : row vector of angular frequencies
%     amps              : row vector of amplitudes (default 1)
%     use_gaussian      : boolean, whether to apply the Gaussian canvas (default true)
%     use_conj          : boolean, whether to keep complex conjugate signal (default false)
%     use_pulse_shaping : boolean, whether to apply pulse shaping (derivative of Gaussian) (default false)
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
    if nargin < 5 || isempty(use_pulse_shaping)
        use_pulse_shaping = false;
    end

    % Allow passing 'pulse_shaping' directly via use_gaussian to avoid updating all wrapper signatures
    if (ischar(use_gaussian) || isstring(use_gaussian)) && strcmpi(use_gaussian, 'pulse_shaping')
        use_pulse_shaping = true;
        use_gaussian = false;
    end

    angular_freqs = angular_freqs(:).';
    amps = amps(:).';

    t_0 = 250;
    T = 100;

    signal = @(t) reshape(conj(sum(amps .* exp(1i * angular_freqs .* t(:)), 2)), size(t));

    if use_pulse_shaping
        signal = @(t) reshape(signal(t) .* exp(-(t(:)-t_0).^2./T^2) .* (-2*(t(:)-t_0)/T^2), size(t));
    elseif use_gaussian
        signal = @(t) reshape(signal(t) .* exp(-(t(:)-t_0).^2./T^2), size(t));
    end
    if ~use_conj
        signal = @(t) real(signal(t));
    end

end
