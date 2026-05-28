function [t_all, x_all, impact_log] = simulate_n_steps( ...
    x0, ...
    params, ...
    nSteps, ...
    controller)

% =========================================================
% SIMULATE_N_STEPS
%
% Simulates multiple hybrid walking steps for the
% RABBIT biped robot.
%
% INPUTS:
%   x0          : initial state
%   params      : robot parameters
%   nSteps      : number of walking steps
%   controller  : optional controller handle
%
% OUTPUTS:
%   t_all       : full simulation time
%   x_all       : full state trajectory
%   impact_log  : impact/reset information
%
% =========================================================

%% --------------------------------------------------------
% Input Handling
% ---------------------------------------------------------

if nargin < 4
    controller = [];
end

fprintf('\n=========================================\n');
fprintf('      MULTI-STEP SIMULATION START\n');
fprintf('=========================================\n');

%% --------------------------------------------------------
% Initialization
% ---------------------------------------------------------

t_all = [];
x_all = [];

impact_log = struct();

x_current = x0;

time_offset = 0;

%% --------------------------------------------------------
% Main Hybrid Simulation Loop
% ---------------------------------------------------------

impact_log(nSteps) = struct();

for step = 1:nSteps
    
    fprintf('\n-----------------------------------------\n');
    fprintf('Simulating Step %d / %d\n', step, nSteps);
    fprintf('-----------------------------------------\n');
    
    %% ----------------------------------------------------
    % Simulate One Continuous Step
    % -----------------------------------------------------
    
    try
        global CURRENT_STEP;
        [t_step, x_step, impact_info] = simulate_one_step( ...
            x_current, ...
            params, ...
            controller);
        CURRENT_STEP = CURRENT_STEP + 1;
        
    catch ME
        
        fprintf('Step simulation failed.\n');
        fprintf('Error:\n%s\n', ME.message);
        
        break;
        
    end
    
    %% ----------------------------------------------------
    % Validate Step
    % -----------------------------------------------------
    
    if isempty(t_step)
        
        fprintf('Empty trajectory returned.\n');
        break;
        
    end
    
    if any(isnan(x_step(:)))
        
        fprintf('NaN detected in simulation.\n');
        break;
        
    end
    
    %% ----------------------------------------------------
    % Shift Time
    % -----------------------------------------------------
    
    t_step = t_step + time_offset;
    
    %% ----------------------------------------------------
    % Concatenate Trajectory
    % -----------------------------------------------------
    
    if isempty(t_all)
        
        t_all = t_step;
        x_all = x_step';
        
    else
        
        % Avoid duplicate transition sample
        t_all = [t_all; t_step(2:end)];
        x_all = [x_all, x_step(2:end,:)'];
        
    end
    
    %% ----------------------------------------------------
    % Store Impact Information
    % -----------------------------------------------------
    
    impact_log(step).step_number = step;
    impact_log(step).impact_time = t_step(end);
    impact_log(step).impact_info = impact_info;
    
    
    %% ----------------------------------------------------
    % Prepare Next Step Initial State
    % -----------------------------------------------------
    
    try
        
        x_current = rabbit_reset_map( ...
            x_step(end,:)', ...
            params);
        
    catch ME
        fprintf('Reset map failed:\n%s\n',ME.message);
        break;
    end
    
    q = x_current(1:7);
    
    [p_st, p_sw, ~, ~, ~, ~] = rabbit_kinematics(q,  packParameters(params));
    
    fprintf('Post-impact stance foot height: %.6f\n', p_st(2));
    fprintf('Post-impact swing foot height : %.6f\n', p_sw(2));
    
    %% ----------------------------------------------------
    % Update Time Offset
    % -----------------------------------------------------
    
    time_offset = t_all(end);
    
    %% ----------------------------------------------------
    % Fall Detection
    % -----------------------------------------------------
    
    nq = round(size(impact_info.state,1)/2);
    q_minus = impact_info.state(1:nq);
    
    [~, p_swing, ~, ~, ~, ~] = rabbit_kinematics(q_minus, packParameters(params));
    
    impact_info.swing_height = p_swing(2);
    
    
    if any(abs(q) > pi)
        
        fprintf('Robot likely fell. Terminating.\n');
        break;
        
    end
    
    fprintf('Step %d completed successfully.\n', step);
    
    fprintf('Step duration: %.3f sec\n', t_step(end)-t_step(1));
    fprintf('Impact swing height: %.6f\n', impact_info.swing_height);
    
end

%% --------------------------------------------------------
% Final Summary
% ---------------------------------------------------------

fprintf('\n=========================================\n');
fprintf('      MULTI-STEP SIMULATION END\n');
fprintf('=========================================\n');

fprintf('Completed Steps : %d\n', length(impact_log));

if ~isempty(t_all)
    
    fprintf('Total Time      : %.3f sec\n', t_all(end));
    fprintf('Trajectory Size : %d samples\n', length(t_all));
    
end

fprintf('=========================================\n');

end
