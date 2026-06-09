function test_normalize_signals
    % TEST_NORMALIZE_SIGNALS Run automated unit tests for normalize_signals.m
    clc;
    fprintf('=== STARTING NORMALIZE_SIGNALS TEST SUITE ===\n\n');
    
    %% Setup Common Test Signals
    % S1: Reference Gaussian peak (Max = 1, Energy approx 1.25)
    S1 = @(t) exp(-t.^2); 
    % S2: Scaled Gaussian peak (Max = 5, Energy is 25x larger)
    S2 = @(t) 5 * exp(-t.^2); 
    % S3: Scaled, complex, shifted Gaussian (Max = 2)
    S3 = @(t) 2i * exp(-(t-1).^2); 
    
    t_test = 0; % Common evaluation point
    
    %% Case 1: Single Function Handle Input (Default Peak)
    fprintf('Test 1: Single function handle input...\n');
    [S2_norm, Rel_N] = normalize_signals(S2);
    
    % Since it is a single signal, its relative factor to "the first signal" (itself) must be 1
    assert(abs(Rel_N - 1) < 1e-6, 'Single handle Rel_N should be 1');
    assert(abs(S2_norm(t_test) - S2(t_test)) < 1e-6, 'Single handle should not scale itself');
    fprintf('  -> Passed.\n\n');

    %% Case 2: Multi-Signal Output Unpacking (Peak Mode)
    fprintf('Test 2: Cell array input with individual variable extraction (Peak)...\n');
    [S1_norm, S2_norm, S3_norm, Rel_N] = normalize_signals({S1, S2, S3}, 'peak');
    
    % S2 has a peak 5x larger than S1. S3 has a peak 2x larger than S1.
    % Rel_N should reflect scaling relative to S1: [1, 5, 2]
    expected_Rel_N = [1, 5, 2];
    assert(max(abs(Rel_N - expected_Rel_N)) < 1e-4, 'Relative peak factors incorrect');
    
    % Once normalized to S1, evaluating them at their respective peaks should yield identical values
    assert(abs(S1_norm(0) - S2_norm(0)) < 1e-4, 'Peak normalization failed between S1 and S2');
    assert(abs(abs(S1_norm(0)) - abs(S3_norm(1))) < 1e-4, 'Peak normalization failed between S1 and S3');
    fprintf('  -> Passed.\n\n');

    %% Case 3: Energy Normalization Mode
    fprintf('Test 3: Energy-based normalization...\n');
    % Setting tight search range around the main bulk of the energy
    [S1_norm_E, S2_norm_E, Rel_N_E] = normalize_signals({S1, S2}, 'energy', [-5, 5]);
    
    % Energy scales with amplitude squared. Amplitude is 5x, so energy is 25x.
    % The RMS/Energy normalization factor ratio should be sqrt(25) = 5.
    assert(max(abs(Rel_N_E - [1, 5])) < 1e-4, 'Relative energy factors incorrect');
    assert(abs(S1_norm_E(0) - S2_norm_E(0)) < 1e-4, 'Energy normalization failed');
    fprintf('  -> Passed.\n\n');

    %% Case 4: Cell Output Fallback (Mismatching nargout)
    fprintf('Test 4: Requesting a single cell array output containing all signals...\n');
    % Asking for 1 output when 3 signals are passed forces the cell wrapper logic
    [norm_cell, Rel_N_cell] = normalize_signals({S1, S2, S3});
    
    assert(iscell(norm_cell), 'Output should fall back to a cell array');
    assert(length(norm_cell) == 3, 'Output cell array missing elements');
    assert(max(abs(Rel_N_cell - [1, 5, 2])) < 1e-4, 'Relative factors inside cell mode incorrect');
    fprintf('  -> Passed.\n\n');

    %% Case 5: Edge Case - Zero/Invalid Signals Error Handling
    fprintf('Test 5: Robustness checks (Dead/Zero signals)...\n');
    S_dead = @(t) zeros(size(t));
    
    % This will trigger the warning blocks in your code; verifying it gracefully handles it
    warning('off', 'all'); % Suppress console spam for expected warnings
    [S1_n, S_dead_n, Rel_N_dead] = normalize_signals({S1, S_dead}, 'peak');
    warning('on', 'all');
    
    % S_dead should fall back to a scale factor of 1 relative to itself, 
    % making its Rel_N factor = 1 / ref_factor = 1 / 1 = 1.
    assert(Rel_N_dead(2) == 1, 'Dead signal fallback factor should be 1');
    assert(S_dead_n(0) == 0, 'Dead signal should remain zero');
    fprintf('  -> Passed.\n\n');

    %% Case 6: Variable Scalar t_guess Expansion
    fprintf('Test 6: Scalar t_guess expansion bounds handling...\n');
    % Test passing a single scalar as t_guess instead of an array
    [~, ~, Rel_N_scalar] = normalize_signals({S1, S2}, 'peak', 0);
    assert(max(abs(Rel_N_scalar - [1, 5])) < 1e-4, 'Scalar t_guess expansion failed');
    fprintf('  -> Passed.\n\n');

    fprintf('=== ALL TESTS PASSED SUCCESSFULLY ===\n');
end