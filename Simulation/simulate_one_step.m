function [t_out, x_out] = simulate_one_step(x0, controller)
% SIMULATE_ONE_STEP Simulates one continuous walking step of the RABBIT robot.

%% Input Handling
arguments
    x0 (:, 1) double {mustBeNonempty, mustBeFinite}
    controller = []
end

%% Solver Configuration
% Increased tolerance and restricted max step for better convergence
options = odeset(...
    'RelTol', 1e-6, ...       
    'AbsTol', 1e-6, ...       
    'MaxStep', 0.01, ...      
    'Events', @rabbit_impact_event);

% Define ODE dynamics
ode_fun = @(t, x) rabbit_ode(t, x, controller);

%% Integration
% Note: Using 0.8 as a hard time limit for a single step
[t_out, x_out, te, xe, ie] = ode45(ode_fun, [0 10], x0, options);

%% Validation
if isempty(t_out) || any(isnan(x_out), 'all')
    error('Simulation failed: Trajectory is empty or contains NaNs.');
end

%% Diagnostics
check_stability(x_out);
end

function check_stability(x_out)
    % Helper to check for physically unrealistic values
    n = size(x_out, 2) / 2;
    if any(abs(x_out(:, 1:n)) > 2*pi, 'all') % Example limit: 2*pi radians
        warning('Unlikely joint angles detected (> 360 degrees).');
    end
    if any(abs(x_out(:, n+1:end)) > 50, 'all') % Example limit: 50 rad/s
        warning('High joint velocities detected.');
    end
end
