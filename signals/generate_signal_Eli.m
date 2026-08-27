function [signal, angular_freqs, amps] = generate_signal_Eli(A, Omega, delta, varargin)
%GENERATE_SIGNAL_ELI Generate the superoscillating signal based on Eli's formula
%
%   [signal, angular_freqs, amps] = generate_signal_Eli(A, Omega, delta, varargin)
%
%   INPUTS:
%     A         : Amplitude parameter (scalar)
%     Omega     : Main angular frequency (scalar)
%     delta     : Angular frequency shift (scalar)
%     varargin  : optional arguments passed directly to generate_signal_base
%                 e.g., use_gaussian, use_conj, use_pulse_shaping
%
%   OUTPUTS:
%     signal        : function handle representing the generated signal
%     angular_freqs : 1x2 vector of frequencies [Omega, Omega - delta]
%     amps          : 1x2 vector of amplitudes [A, -(A-1)]
%
%   FORMULA:
%     Signal(t) = A * exp(i * Omega * t) - (A - 1) * exp(i * (Omega - delta) * t)

    % Default values if not provided
    if nargin < 1 || isempty(A)
        A = 10;
    end
    if nargin < 2 || isempty(Omega)
        Omega = 1;
    end
    if nargin < 3 || isempty(delta)
        delta = 0.1;
    end

    % Construct the frequencies and amplitudes according to the formula
    amps = [A, -(A - 1)];
    angular_freqs = [Omega, Omega - delta];

    % Call the base generator
    [signal, angular_freqs, amps] = generate_signal_base(angular_freqs, amps, varargin{:});

end
