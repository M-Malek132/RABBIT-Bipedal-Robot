function [t_all, x_all, impact_log] = simulate_n_steps(x0, nSteps, controller)
% SIMULATE_N_STEPS Simulates multiple hybrid walking steps for RABBIT.

%% Input Handling
arguments
    x0 (:, 1) double {mustBeNonempty, mustBeFinite}
    nSteps (1, 1) double {mustBeInteger, mustBePositive}
    controller = []
end

fprintf('\n=== Multi-Step Simulation Start ===\n');

%% Initialization
t_all = [];
x_all = [];
impact_log = struct('step_number', {}, 'impact_time', {}, 'impact_info', {}, ...
                    'pre_impact', {}, 'post_impact', {});

x_current = x0;
time_offset = 0;

%% Simulation Loop
for step = 1:nSteps
    fprintf('Simulating Step %d / %d...\n', step, nSteps);

    % Simulate
    try
        [t_step, x_step, impact_info] = simulate_one_step(x_current, controller);
    catch ME
        fprintf('Error at step %d: %s\n', step, ME.message);
        break;
    end

    % Validation
    if isempty(t_step) || any(isnan(x_step), 'all')
        error('Simulation returned invalid trajectory at step %d.', step);
    end

    % Concatenate (Offset time to maintain continuous timeline)
    if isempty(t_all)
        t_all = t_step;
        x_all = x_step;
    else
        t_all = [t_all; t_step(2:end) + time_offset];
        x_all = [x_all; x_step(2:end, :)];
    end

    % Store Impact Data
    impact_log(step).step_number = step;
    impact_log(step).impact_time = t_all(end);
    impact_log(step).impact_info = impact_info;
    impact_log(step).pre_impact  = x_step(end, :)';

    % Reset Map
    x_current = rabbit_reset_map(impact_log(step).pre_impact);
    impact_log(step).post_impact = x_current;

    % Update Tracking
    time_offset = t_all(end);
    fprintf('Step %d success. Duration: %.3fs\n', step, t_step(end) - t_step(1));
end

%% Final Summary
fprintf('=== Simulation End: Completed %d steps ===\n', numel(impact_log));
end
