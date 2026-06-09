% Add all project subfolders to search path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

test_normalize_signals();
test_compute_instantaneous_frequency();