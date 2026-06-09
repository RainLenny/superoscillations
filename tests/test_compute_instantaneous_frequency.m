function test_compute_instantaneous_frequency
    % TEST_COMPUTE_INSTANTANEOUS_FREQUENCY Run automated unit tests for compute_instantaneous_frequency.m
    clc;
    fprintf('=== STARTING COMPUTE_INSTANTANEOUS_FREQUENCY TEST SUITE ===\n\n');
    
    %% Setup path
    % Ensure functions are in path if they aren't already
    addpath(genpath(fileparts(mfilename('fullpath'))));
    
    %% Test 1: Pure cosine wave with known frequency
    fprintf('Test 1: Pure cosine wave with known frequency...\n');
    omega_c = 2.5;
    t = (0:0.001:10)';
    
    inst_freq = compute_instantaneous_frequency(omega_c, 1, t);
    
    % The analytical instantaneous frequency of a pure cosine should be exactly omega_c everywhere
    
    mean_freq = mean(inst_freq);
    max_dev = max(abs(inst_freq - omega_c));
    
    fprintf('  Expected freq: %.4f, Mean calculated freq: %.4f, Max deviation: %.6e\n', ...
        omega_c, mean_freq, max_dev);
    
    assert(abs(mean_freq - omega_c) < 1e-10, 'Mean instantaneous frequency dev too large');
    assert(max_dev < 1e-10, 'Max deviation from expected carrier frequency too large');
    fprintf('  -> Passed.\n\n');
    
    %% Test 2: Input with complex amplitudes
    fprintf('Test 2: Complex amplitudes (phase offset)...\n');
    amp_complex = 1 * exp(1i * pi/4);
    
    inst_freq_complex = compute_instantaneous_frequency(omega_c, amp_complex, t);
    assert(max(abs(inst_freq_complex - omega_c)) < 1e-10, ...
        'Complex amplitude did not correctly yield carrier frequency');
    fprintf('  -> Passed.\n\n');
    
    fprintf('=== ALL TESTS PASSED SUCCESSFULLY ===\n');
end
