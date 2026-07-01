function [t_out, x_out, impact_info] = simulate_one_step(x0, controller)
% SIMULATE_ONE_STEP Simulates one continuous walking step of the RABBIT robot.
%
%   [t_out, x_out, impact_info] = simulate_one_step(x0, controller)
%
%   INPUTS:
%       x0          - Initial state [q; dq]
%       controller  - (Optional) Controller function handle
%
%   OUTPUTS:
%       t_out       - Simulation time vector
%       x_out       - State trajectory
%       impact_info - Struct containing event detection details

    % Input validation
    narginchk(1, 2);
    if nargin < 2, controller = []; end

    fprintf('Starting single-step simulation...\n');

    %% Solver Configuration
    % Note: rabbit_impact_event and rabbit_ode must now access 
    % parameters via global variables or internal constants.
    options = odeset(...
        'RelTol', 1e-3, ...       
        'AbsTol', 1e-3, ...       
        'MaxStep', 0.05, ...      
        'Events', @(t,x) rabbit_impact_event(t, x));

    ode_fun = @(t,x) rabbit_ode(t, x, controller);

    %% Integrate Dynamics
    try
        [t_out, x_out, te, xe, ie] = ode45(ode_fun, [0 0.8], x0, options);
    catch ME
        % Safely handle the global step counter
        if exist('CURRENT_STEP', 'var')
            fprintf(2, 'Integration failed at step %d: %s\n', CURRENT_STEP, ME.message);
        else
            fprintf(2, 'Integration failed: %s\n', ME.message);
        end
        rethrow(ME);
    end

    %% Post-Integration Validation
    if isempty(t_out) || any(isnan(x_out(:)))
        error('Simulation failed: Trajectory is empty or contains NaNs.');
    end

    %% Process Impact Info
    impact_info = struct('detected', false, 'time', [], 'state', [], 'index', []);
    
    if ~isempty(te)
        fprintf('Impact detected at t = %.4f sec\n', te(end));
        impact_info.detected = true;
        impact_info.time     = te(end);
        impact_info.state    = xe(end,:)';
        impact_info.index    = ie(end);
    else
        fprintf('No impact detected.\n');
    end

    %% Diagnostics
    fprintf('Simulation finished in %.4f sec (%d points)\n', t_out(end), length(t_out));
    check_stability(x_out);
end

function check_stability(x_out)
    % Helper to check for physically unrealistic values
    n = size(x_out, 2) / 2;
    if any(abs(x_out(:, 1:n)) > 10, 'all')
        warning('Large joint angles detected.');
    end
    if any(abs(x_out(:, n+1:end)) > 100, 'all')
        warning('Large joint velocities detected.');
    end
end
