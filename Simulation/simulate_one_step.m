function [t_out, x_out, impact_info] = simulate_one_step( ...
    x0, ...
    params, ...
    controller)

% =========================================================
% SIMULATE_ONE_STEP
%
% Simulates one continuous walking step of the
% RABBIT robot until impact occurs.
%
% INPUTS:
%   x0          : initial state [q; dq]
%   params      : robot parameter structure
%   controller  : optional controller function handle
%
% OUTPUTS:
%   t_out       : simulation time vector
%   x_out       : state trajectory
%   impact_info : impact event information
%
% =========================================================

%% --------------------------------------------------------
% Input Handling
% ---------------------------------------------------------

if nargin < 3
    controller = [];
end

fprintf('\nStarting single-step simulation...\n');

%% --------------------------------------------------------
% Simulation Settings
% ---------------------------------------------------------

tspan = [0 5];

options = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'Events', @(t,x) rabbit_impact_event(t, x, params));

%% --------------------------------------------------------
% ODE Function
% ---------------------------------------------------------

ode_fun = @(t,x) rabbit_ode( ...
    t, ...
    x, ...
    params, ...
    controller);

%% --------------------------------------------------------
% Integrate Dynamics
% ---------------------------------------------------------

try

    [t_out, x_out, te, xe, ie] = ode45( ...
        ode_fun, ...
        tspan, ...
        x0, ...
        options);

catch ME

    fprintf('ODE integration failed.\n');
    rethrow(ME);

end

%% --------------------------------------------------------
% Validate Output
% ---------------------------------------------------------

if isempty(t_out)
    error('Empty trajectory returned from ODE solver.');
end

if any(isnan(x_out(:)))
    error('NaN detected during integration.');
end

%% --------------------------------------------------------
% Impact Information
% ---------------------------------------------------------

impact_info = struct();

if isempty(te)

    fprintf('No impact detected.\n');

    impact_info.detected = false;
    impact_info.time     = [];
    impact_info.state    = [];
    impact_info.index    = [];

else

    fprintf('Impact detected at t = %.4f sec\n', te(end));

    impact_info.detected = true;
    impact_info.time     = te(end);
    impact_info.state    = xe(end,:)';
    impact_info.index    = ie(end);

end

%% --------------------------------------------------------
% Trajectory Diagnostics
% ---------------------------------------------------------

fprintf('Simulation duration : %.4f sec\n', t_out(end));
fprintf('Trajectory samples  : %d\n', length(t_out));

%% --------------------------------------------------------
% Optional Stability Check
% ---------------------------------------------------------

n = size(x_out,2)/2;

q  = x_out(:,1:n);
dq = x_out(:,n+1:end);

if any(abs(q(:)) > 10)

    warning('Large joint angles detected.');

end

if any(abs(dq(:)) > 100)

    warning('Large joint velocities detected.');

end

fprintf('Single-step simulation completed.\n');

end


