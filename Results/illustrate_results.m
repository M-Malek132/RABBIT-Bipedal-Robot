function illustrate_results(result_file)
% ILLUSTRATE_RESULTS Load a simulation result MAT file and visualize it.
%
% Usage:
%   illustrate_results('result/simulation_results_YYYY-mm-dd_HH-MM-SS.mat')
%
% The function:
%   1) loads t_all, x_all, params, impact_log
%   2) plots joint trajectories
%   3) plots base position and velocity
%   4) animates the robot if animation function exists
%
% Inputs:
%   result_file : path to .mat file

if nargin < 1 || isempty(result_file)
    error('Please provide the path to a saved result .mat file.');
end

if ~exist(result_file, 'file')
    error('Result file not found: %s', result_file);
end

S = load(result_file);

% Basic checks
if ~isfield(S, 't_all') || ~isfield(S, 'x_all')
    error('MAT file must contain t_all and x_all.');
end

t_all = S.t_all;
x_all = (S.x_all)';

if isfield(S, 'params')
    params = S.params;
else
    params = [];
end

fprintf('Loaded results from: %s\n', result_file);

% -----------------------------
% Plot state trajectories
% -----------------------------
nq = round(size(x_all, 2) / 2);
if mod(size(x_all,2),2) ~= 0
    error('x_all must have an even number of columns (q and dq).');
end

q_all  = x_all(:, 1:nq);
dq_all = x_all(:, nq+1:end);

figure('Name','Rabbit Results','Color','w');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

% Base position x, z, theta
nexttile;
plot(t_all, q_all(:,1), 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Base x');
title('Base Horizontal Position');

nexttile;
plot(t_all, q_all(:,2), 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Base z');
title('Base Vertical Position');

nexttile;
plot(t_all, q_all(:,3), 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('Torso angle');
title('Torso Angle');

% Joint angles
nexttile;
plot(t_all, q_all(:,4:end), 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Angle [rad]');
title('Joint Angles');
legend({'q1','q2','q3','q4'}, 'Location','best');

% Joint velocities
nexttile;
plot(t_all, dq_all(:,4:end), 'LineWidth', 1.2); grid on;
xlabel('Time [s]'); ylabel('Velocity [rad/s]');
title('Joint Velocities');
legend({'dq1','dq2','dq3','dq4'}, 'Location','best');

% Optional: total speed norm
nexttile;
plot(t_all, vecnorm(dq_all,2,2), 'LineWidth', 1.5); grid on;
xlabel('Time [s]'); ylabel('||dq||');
title('Velocity Norm');

% -----------------------------
% Impact log (if available)
% -----------------------------
if isfield(S, 'impact_log')
    figure('Name','Impact Log','Color','w');
    if isnumeric(S.impact_log) && ~isempty(S.impact_log)
        stem(S.impact_log, ones(size(S.impact_log)), 'filled');
        xlabel('Impact index / time entry');
        ylabel('Impact event');
        title('Impact Log');
        grid on;
    else
        text(0.1,0.5,'impact_log exists but is not numeric / is empty');
        axis off;
    end
end

% -----------------------------
% Animate if possible
% -----------------------------
if ~isempty(params)
    try
        figure('Name','Rabbit Animation','Color','w');
        animate_rabbit_stepping_stones(x_all, params);
    catch ME
        warning('Animation failed: %s', ME.message);
    end
else
    fprintf('params not found in MAT file; skipping animation.\n');
end

% Add this to illustrate_results.m to visualize when your constraints break
if isfield(S, 'delta_log') % Assuming you save slacks
    figure('Name','Constraint Slack Variables');
    plot(t_all, S.delta_log(:,1), 'r', 'LineWidth', 1.5); hold on;
    plot(t_all, S.delta_log(:,2), 'b', 'LineWidth', 1.5);
    legend('Delta CLF', 'Delta CBF');
    title('Constraint Slack (Violation Magnitude)');
    grid on;
end

end
