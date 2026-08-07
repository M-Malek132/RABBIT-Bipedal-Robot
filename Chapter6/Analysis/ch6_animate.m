function fig = ch6_animate(sim, p, file, n_frames)
%CH6_ANIMATE  Stick-figure snapshots over the stepping stones.  Figs. 6.6, 6.27.
%
%   fig = ch6_animate(sim, p)
%   fig = ch6_animate(sim, p, file, n_frames)
%
% Draws the robot at n_frames poses per step, overlaid on the terrain, with the
% stones drawn as the black pillars of Fig. 6.6 -- each pillar spanning the
% window [l_min, l_max] the step was actually asked to land in, placed at the
% world position of the stance foot it was measured from.
%
% Placing the pillars in WORLD coordinates is the only part that needs care.
% l_min and l_max are measured from the stance foot of their own step, and the
% stance foot advances every step, so the world position of stone k is
%
%       x_stance(k) + [l_min(k), l_max(k)]
%
% where x_stance(k) is read from P_st at the start of step k. Drawing them at a
% fixed spacing instead would produce a figure in which the robot appears to hit
% stones it missed.
%
% Inputs
%   sim      : output of ch6_simulate
%   p        : parameter struct
%   file     : optional PNG path
%   n_frames : poses drawn per step (default 6)
%
% See also CH6_SIMULATE, CH3_BODY_POINTS, CH6_PLOT_STONES.

if nargin < 3, file = ''; end
if nargin < 4 || isempty(n_frames), n_frames = 6; end

n = sim.n_ok;
if n == 0
    error('ch6_animate:empty', 'The run completed no steps.');
end

fig = figure('Color', 'w', 'Position', [60 200 1300 420]);
hold on; box on; axis equal;

fade = @(a) [0.30 0.30 0.35] * a + [1 1 1] * (1 - a);

for k = 1:n
    s  = sim.steps(k);
    q0 = s.x(1:p.nq, 1);
    x_st = P_st(q0);

    % ---- the stone: a pillar under the admissible landing window
    lo = x_st(1) + sim.l_min(k);
    hi = x_st(1) + sim.l_max(k);
    patch([lo hi hi lo], [-0.10 -0.10 0 0], [0.15 0.15 0.15], ...
          'EdgeColor', 'none');

    % ---- the robot, faded from light (start of step) to dark (impact)
    idx = round(linspace(1, size(s.x, 2), n_frames));
    for f = 1:numel(idx)
        q = s.x(1:p.nq, idx(f));
        draw_robot(q, fade(0.25 + 0.75*f/numel(idx)));
    end

    % ---- where the foot actually landed
    q_end = s.x(1:p.nq, end);
    fe = P_sw(q_end);
    plot(fe(1), fe(2), 'o', 'MarkerSize', 7, 'MarkerEdgeColor', 'k', ...
         'MarkerFaceColor', hit_color(sim.hit(k)));
end

yline(0, 'k-', 'LineWidth', 1.2);
xlabel('x (m)'); ylabel('z (m)');
title(sprintf(['RABBIT over %d discrete footholds  --  %d/%d placements on ' ...
               'the stone'], n, sum(sim.hit), n));
ylim([-0.12, 1.35]);
grid on;

if ~isempty(file)
    exportgraphics(fig, file, 'Resolution', 150);
    fprintf('[ch6_animate] wrote %s\n', file);
end

end

% ---------------------------------------------------------------------------
function draw_robot(q, col)
pts = ch3_body_points(q);
seg = {[pts.hip, pts.torso_top], ...
       [pts.hip, pts.stance_knee], [pts.stance_knee, pts.stance_foot], ...
       [pts.hip, pts.swing_knee],  [pts.swing_knee,  pts.swing_foot]};
w   = [2.6, 1.8, 1.8, 1.8, 1.8];
for i = 1:numel(seg)
    plot(seg{i}(1,:), seg{i}(2,:), '-', 'Color', col, 'LineWidth', w(i));
end
end

function c = hit_color(b)
if b, c = [0.15 0.65 0.20]; else, c = [0.85 0.15 0.15]; end
end
