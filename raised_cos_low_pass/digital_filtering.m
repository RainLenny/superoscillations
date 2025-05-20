generate_signals

f_sampling = 60/(2*pi);      % sampling rate in Hz
window_length_for_filter = inf;

digital_filtering_func(sinc_signal,superoscillations_signal,angular_freqs,...
    filter_signal, f_sampling,window_length_for_filter)