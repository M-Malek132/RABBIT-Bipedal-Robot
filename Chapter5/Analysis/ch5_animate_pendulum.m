function gif = ch5_animate_pendulum(runs, gifpath, varargin)
%CH5_ANIMATE_PENDULUM  Side-by-side animation of the Fig. 5.4 runs.
%
%   gif = ch5_animate_pendulum(runs, gifpath)
%   gif = ch5_animate_pendulum(runs, gifpath, 'fps', 25, 'stride', 40)
%
% One panel per run, all stepping in lockstep, with the p2min line drawn in
% each. The tables say the baseline "violated by 0.98 m"; this is where you see
% what that means -- the unconstrained arm sweeps straight down through the
% line while the constrained ones fold and slide along it.
%
% Watching the ELBOW is the point. The end-effector traces look similar until
% the constraint bites; the difference is in theta2, and that is a joint angle
% rather than a plotted signal.
%
% Inputs
%   runs    : struct array with .sim and .label
%   gifpath : output .gif path
%   'fps'    frames per second in the GIF (default 25)
%   'stride' take every Nth control sample (default: aim for ~fps*T frames)
%
% Output
%   gif : the path written, or '' if nothing was
%
% See also CH5_PLOT_PENDULUM, CH5_PEND_POINTS, CH5_MAIN.

o = struct('fps', 25, 'stride', []);
for k = 1:2:numel(varargin), o.(lower(varargin{k})) = varargin{k+1}; end

n = numel(runs);
if n == 0, gif = ''; return; end

s1 = runs(1).sim;
T  = s1.t(end);
if isempty(o.stride)
    o.stride = max(1, round(numel(s1.t) / (o.fps * T)));
end
idx = 1:o.stride:numel(s1.t);

reach = sum(s1.p.plant.l);

fig = figure('Name', 'ch5: pendulum, side by side', ...
             'Position', [60 60 380*n 440], 'Color', 'w');
ax = gobjects(1, n);
for i = 1:n
    ax(i) = subplot(1, n, i);
    hold(ax(i), 'on'); grid(ax(i), 'on'); axis(ax(i), 'equal');
    xlim(ax(i), [-1.2 1.2]*reach); ylim(ax(i), [-1.2 1.2]*reach);
    yline(ax(i), runs(i).sim.p.constraint.value, 'r-', 'LineWidth', 1.5);
    title(ax(i), runs(i).label, 'Interpreter', 'none', 'FontWeight', 'normal');
end

gif = gifpath;
first = true;

for k = idx
    for i = 1:n
        s = runs(i).sim;
        kk = min(k, numel(s.t));

        delete(findobj(ax(i), 'Tag', 'arm'));

        P = ch5_pend_points(s.x(1:2, kk), s.p);
        plot(ax(i), squeeze(P(1,:,1)), squeeze(P(2,:,1)), '-o', ...
             'Color', [0.1 0.1 0.1], 'LineWidth', 3, 'MarkerSize', 6, ...
             'MarkerFaceColor', [0.1 0.1 0.1], 'Tag', 'arm');

        % colour the tip by whether the barrier is binding right now
        if s.cbf_active(kk), c = [0.85 0.1 0.1]; else, c = [0 0.5 0.1]; end
        plot(ax(i), P(1,3,1), P(2,3,1), 'o', 'MarkerSize', 10, ...
             'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'Tag', 'arm');

        xlabel(ax(i), sprintf('t = %5.2f s    p^y_2 = %+.3f m', ...
                              s.t(kk), s.h(kk) + s.p.constraint.value));
    end

    drawnow limitrate;
    frame = getframe(fig);
    [A, map] = rgb2ind(frame2im(frame), 256);

    if first
        imwrite(A, map, gif, 'gif', 'LoopCount', inf, 'DelayTime', 1/o.fps);
        first = false;
    else
        imwrite(A, map, gif, 'gif', 'WriteMode', 'append', 'DelayTime', 1/o.fps);
    end
end

close(fig);

end
