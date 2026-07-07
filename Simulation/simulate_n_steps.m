function [t_all, x_all] = simulate_n_steps(x0, nSteps, controller)
% SIMULATE_N_STEPS Simulates multiple hybrid walking steps for RABBIT.

fprintf('\n=== Multi-Step Simulation Start ===\n');

t_all = [];
x_all = [];

x_current = x0;
time_offset = 0;

for step = 1:nSteps
    fprintf('Simulating Step %d / %d...\n', step, nSteps);

    % 1. Simulate the continuous phase
    % simulate_one_step should return the trajectory up to the point of impact
    [t_step, x_step] = simulate_one_step(x_current, controller);

    if isempty(t_step) || any(isnan(x_step), 'all')
        fprintf('Simulation failed or terminated at step %d.\n', step);
        break;
    end

    % 2. Process the hybrid transition
    x_minus = x_step(end, :)';
    
    % A. Physical Impact: Calculate velocity jump (Momentum balance)
    x_after_impact = rabbit_impact_map(x_minus);
    
    % B. Reset/Relabel: Swap legs and shift coordinates (Reset map logic)
    x_next_start = rabbit_reset_map(x_after_impact);

    % 3. Store Data
    if isempty(t_all)
        t_all = t_step;
        x_all = x_step;
    else
        t_all = [t_all; t_step(2:end) + time_offset];
        x_all = [x_all; x_step(2:end, :)];
    end

    % 4. Prepare for next iteration
    x_current = x_next_start;
    time_offset = t_all(end);
    
    fprintf('Step %d success. Duration: %.3fs\n', step, t_step(end) - t_step(1));
end
end
