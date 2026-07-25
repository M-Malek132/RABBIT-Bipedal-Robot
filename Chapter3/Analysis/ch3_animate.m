function ch3_animate(x0, alpha, p, n_steps, gifpath)
%CH3_ANIMATE  Animate the gait and optionally write a GIF.
%
%   ch3_animate(x0, alpha, p, n_steps)
%   ch3_animate(x0, alpha, p, n_steps, 'Results/ch3_walk.gif')
%
% Simulates n_steps under p.controller and draws the stick figure, stance leg
% solid and swing leg dashed, with the ground line and the accumulated
% footprints.
%
% This is simulation, not the collocation trajectory -- so if the gait was
% built on too coarse a mesh, this is where it visibly falls over. Worth
% watching for that reason alone.
%
% Requires a normal (JVM-enabled) MATLAB session.
%
% See also CH3_SIMULATE, CH3_PLOT_GAIT, CH3_BODY_POINTS.

if nargin < 4 || isempty(n_steps), n_steps = 5; end
if nargin < 5, gifpath = ''; end

sim = ch3_simulate(x0, alpha, p, n_steps);
if sim.n_ok == 0
    error('ch3_animate:noSteps', 'No steps completed: %s', sim.reason);
end
fprintf('ch3_animate: %d/%d steps (%s)\n', sim.n_ok, n_steps, sim.reason);

X = sim.x;
nt = size(X, 2);
stride = max(1, round(nt / 200));      % ~200 frames regardless of solver grid
frames = 1:stride:nt;

xs = X(1,:);
fig = figure('Color','w','Position',[100 100 900 420]);
ax  = axes(fig); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlim(ax, [min(xs)-1.0, max(xs)+1.0]);
ylim(ax, [-0.15, 1.6]);
xlabel(ax,'x [m]'); ylabel(ax,'z [m]');

first = true;
for k = frames
    cla(ax);
    plot(ax, xlim(ax), [0 0], 'k-', 'LineWidth', 2);

    b = ch3_body_points(X(1:p.nq, k));
    plot(ax, [b.hip(1) b.torso_top(1)], [b.hip(2) b.torso_top(2)], 'k-', 'LineWidth', 3);
    plot(ax, [b.hip(1) b.stance_knee(1) b.stance_foot(1)], ...
             [b.hip(2) b.stance_knee(2) b.stance_foot(2)], '-o', ...
             'Color',[0.85 0.2 0.2], 'LineWidth', 2.5, 'MarkerSize', 5);
    plot(ax, [b.hip(1) b.swing_knee(1) b.swing_foot(1)], ...
             [b.hip(2) b.swing_knee(2) b.swing_foot(2)], '--s', ...
             'Color',[0.2 0.4 0.85], 'LineWidth', 2, 'MarkerSize', 5);
    plot(ax, b.hip(1), b.hip(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);

    title(ax, sprintf('t = %.2f s   |   %s   |   %d steps', ...
                      sim.t(k), p.controller, sim.n_ok));
    drawnow limitrate;

    if ~isempty(gifpath)
        fr = getframe(fig);
        [ind, cm] = rgb2ind(frame2im(fr), 256);
        if first
            imwrite(ind, cm, gifpath, 'gif', 'LoopCount', Inf, 'DelayTime', 0.04);
            first = false;
        else
            imwrite(ind, cm, gifpath, 'gif', 'WriteMode', 'append', 'DelayTime', 0.04);
        end
    end
end

if ~isempty(gifpath)
    fprintf('Wrote %s\n', gifpath);
end

end
