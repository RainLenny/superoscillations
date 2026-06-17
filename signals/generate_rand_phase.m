function [Rand_signal, angular_freqs_SO] = generate_rand_phase(input_freqs, input_amps, use_gaussian, use_conj, random_seed)
%GENERATE_RAND_PHASE Generate a signal with given amplitudes/frequencies but random phase
%
%   [Rand_signal, angular_freqs_SO] = generate_rand_phase(input_freqs, input_amps, use_gaussian, use_conj, random_seed)
%
%   INPUTS:
%     input_freqs   : vector of angular frequencies
%     input_amps    : vector of corresponding amplitudes
%     use_gaussian  : boolean, whether to apply the Gaussian canvas (default true)
%     use_conj      : boolean, whether to keep complex conjugate signal (default false)
%     random_seed   : integer, seed for the random number generator (default 1)
%
%   OUTPUTS:
%     Rand_signal       : function handle of the randomized phase signal
%     angular_freqs_SO  : 1xN vector of angular frequencies used

    if nargin < 3 || isempty(use_gaussian)
        use_gaussian = true;
    end
    if nargin < 4 || isempty(use_conj)
        use_conj = false;
    end
    if nargin < 5 || isempty(random_seed)
        random_seed = 1;
    end

    angular_freqs_SO = input_freqs;
    N_SO = length(angular_freqs_SO);
    
    rng(random_seed);
    tau_rand = rand(1, N_SO);
    amps_rand = abs(input_amps) .* exp(-1i .* tau_rand);
    
    [Rand_signal, angular_freqs_SO] = generate_signal_base(angular_freqs_SO, amps_rand, use_gaussian, use_conj);

end
