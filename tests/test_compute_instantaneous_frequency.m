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
    y = cos(omega_c * t);
    
    [inst_freq, phase, z] = compute_instantaneous_frequency(y, t);
    
    % The derivative of the phase of cos(w*t) analytic signal should be exactly w.
    % We exclude the boundary elements because gradient uses one-sided difference at edges.
    mid_indices = (t > 1) & (t < 9);
    
    mean_freq = mean(inst_freq(mid_indices));
    max_dev = max(abs(inst_freq(mid_indices) - omega_c));
    
    fprintf('  Expected freq: %.4f, Mean calculated freq: %.4f, Max deviation in middle: %.6e\n', ...
        omega_c, mean_freq, max_dev);
    
    assert(abs(mean_freq - omega_c) < 5e-3, 'Mean instantaneous frequency dev too large');
    assert(max_dev < 3e-2, 'Max deviation from expected carrier frequency too large');
    fprintf('  -> Passed.\n\n');
    
    %% Test 2: Input with complex values (should operate on real part)
    fprintf('Test 2: Complex input (should yield same result as real part)...\n');
    y_complex = y + 1i * sin(omega_c * t);
    
    inst_freq_complex = compute_instantaneous_frequency(y_complex, t);
    assert(max(abs(inst_freq_complex(mid_indices) - omega_c)) < 3e-2, ...
        'Complex input did not correctly extract real part and calculate frequency');
    fprintf('  -> Passed.\n\n');
    
    fprintf('=== ALL TESTS PASSED SUCCESSFULLY ===\n');
end
