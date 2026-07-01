%% SETUP & INITIALIZATION
clc; clear; close all;
fprintf("Testing simulate_one_step...\n");

startup; % Initialize project paths
% One-time initialization at the start of your program
config('init');

% Then proceed to simulation without passing params
x0 = make_initial_state(); % Ensure this also uses config('get') internally
simulate_one_step(x0, @rabbit_controller);

%% SIMULATION
% We define the controller as a closure to avoid passing 'params' repeatedly
controller = @(t, x) rabbit_controller(t, x, params);

% If you truly want to remove 'params' from the signature of simulate_one_step:
% Ensure your global or persistent settings are updated here.
[t, x, impact] = simulate_one_step(x0, controller);

%% VISUALIZATION & DIAGNOSTICS
% Plotting
figure; plot(t, x(:, 1:nq), 'LineWidth', 1.5); title("Joint Angles"); grid on;
figure; plot(t, x(:, nq+1:end), 'LineWidth', 1.5); title("Joint Velocities"); grid on;

% Impact Data
disp("Impact information:"); disp(impact);

% Torque Calculation
u = arrayfun(@(k) controller(t(k), x(k, :)'), 1:length(t), 'UniformOutput', false);
u = cell2mat(u)';
fprintf("Maximum controller torque: %f Nm\n", max(abs(u(:))));

figure; plot(t, u, 'LineWidth', 1.5); title("Controller Effort"); 
legend("St. Hip", "St. Knee", "Sw. Hip", "Sw. Knee"); grid on;

% Animation
animate_rabbit_stepping_stones(x', params);

% Final Verification
[stance_final, ~] = rabbit_kinematics(x(end, 1:7), params_packed);
fprintf("Final stance foot position:\n"); disp(stance_final);
