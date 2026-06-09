function [SO_signal, angular_freqs_SO] = generate_SO_from_dat_file(dat_file_name, freq_scaling, amp_scaling, use_gaussian, use_conj)
%GENERATE_SO_FROM_DAT_FILE Generate the Superoscillating (SO) signal from a .mat data file.
%
%   [SO_signal, angular_freqs_SO] = generate_SO_from_dat_file(dat_file_name, freq_scaling, amp_scaling, use_gaussian, use_conj)
%
%   INPUTS:
%     dat_file_name : name of the data file in the signals/data folder (e.g. 'SO_Baranov' or 'SO_max_0dot7')
%                     If empty or omitted, defaults to 'SO_Baranov.mat'.
%                     Automatically appends '.mat' extension if not present.
%                     Performs case-insensitive matching.
%     freq_scaling  : multiplier for frequencies (default 1)
%     amp_scaling   : multiplier for amplitudes (default 1)
%     use_gaussian  : boolean, whether to apply the Gaussian canvas (default true)
%     use_conj      : boolean, whether to keep complex conjugate signal (default false)
%
%   OUTPUTS:
%     SO_signal         : function handle of the SO signal
%     angular_freqs_SO  : 1xN vector of scaled angular frequencies

    % Default filename if not provided or empty
    if nargin < 1 || isempty(dat_file_name)
        dat_file_name = 'SO_Baranov.mat';
    end

    % Ensure it is a character array
    dat_file_name = char(dat_file_name);

    % Append .mat if not present
    if ~endsWith(dat_file_name, '.mat', 'IgnoreCase', true)
        dat_file_name = [dat_file_name, '.mat'];
    end

    % Find the absolute path to the data folder
    current_dir = fileparts(mfilename('fullpath'));
    data_dir = fullfile(current_dir, 'data');

    % Search directory case-insensitively for the matching file
    files = dir(data_dir);
    matched_filename = '';
    for k = 1:length(files)
        if ~files(k).isdir && strcmpi(files(k).name, dat_file_name)
            matched_filename = files(k).name;
            break;
        end
    end

    if isempty(matched_filename)
        error('generate_SO_from_dat_file:FileNotFound', 'Signal data file "%s" not found in folder "%s".', dat_file_name, data_dir);
    end

    % Load the mat file
    full_path = fullfile(data_dir, matched_filename);
    data = load(full_path, 'amps_SO', 'angular_freqs_SO');

    if ~isfield(data, 'amps_SO') || ~isfield(data, 'angular_freqs_SO')
        error('generate_SO_from_dat_file:InvalidData', 'Loaded file "%s" is missing "amps_SO" or "angular_freqs_SO".', matched_filename);
    end

    if nargin < 2 || isempty(freq_scaling)
        freq_scaling = 1;
    end
    if nargin < 3 || isempty(amp_scaling)
        amp_scaling = 1;
    end
    if nargin < 4 || isempty(use_gaussian)
        use_gaussian = true;
    end
    if nargin < 5 || isempty(use_conj)
        use_conj = false;
    end

    amps_SO = data.amps_SO * amp_scaling;
    angular_freqs_SO = data.angular_freqs_SO * freq_scaling;

    [SO_signal, angular_freqs_SO] = generate_signal_base(angular_freqs_SO, amps_SO, use_gaussian, use_conj);

end
