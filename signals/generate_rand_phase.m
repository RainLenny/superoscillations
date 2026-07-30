function [Rand_signal, angular_freqs_SO, amps_rand] = generate_rand_phase(input_freqs, input_amps, varargin)
%GENERATE_RAND_PHASE Generate a signal with given amplitudes/frequencies but random phase
%
%   [Rand_signal, angular_freqs_SO] = generate_rand_phase(input_freqs, input_amps, varargin)
%
%   INPUTS:
%     input_freqs   : vector of angular frequencies
%     input_amps    : vector of corresponding amplitudes
%     varargin      : optional arguments:
%                       {1} use_gaussian (default true)
%                       {2} use_conj (default false)
%                       {3} random_seed (default 1)
%                       {4+} use_pulse_shaping, etc. (passed to base)
%
%   OUTPUTS:
%     Rand_signal       : function handle of the randomized phase signal
%     angular_freqs_SO  : 1xN vector of angular frequencies used

    random_seed = 1;
    if length(varargin) >= 3 && ~isempty(varargin{3})
        random_seed = varargin{3};
    end
    
    % Prepare arguments for base generator by removing random_seed
    base_varargin = varargin;
    if length(base_varargin) >= 3
        base_varargin(3) = [];
    end

    angular_freqs_SO = input_freqs;
    N_SO = length(angular_freqs_SO);
    
    rng(random_seed);
    tau_rand = 2*pi*rand(1, N_SO);
    amps_rand = abs(input_amps) .* exp(-1i .* tau_rand);
    
    [Rand_signal, angular_freqs_SO , amps_rand] = generate_signal_base(angular_freqs_SO, amps_rand, base_varargin{:});

end
