function b = ch4_draw_robot(ax, q, dimmed)
%CH4_DRAW_ROBOT  The stick figure, in ch3_animate's colours.
%
%   b = ch4_draw_robot(ax, q)
%   b = ch4_draw_robot(ax, q, dimmed)
%
% Shared by the race animation (ch4_animate) and the load-study snapshots
% (ch4_plot_load), so the robot looks the same in both.
%
% A finished run is drawn washed out rather than removed: the pose still carries
% information (where it ended up, and in what attitude), but it must not read as
% a robot that is still walking.
%
% Output
%   b : the body points drawn (ch3_body_points), so a caller can place things
%       on the robot -- a carried load at b.hip, for one
%
% See also CH4_ANIMATE, CH4_PLOT_LOAD, CH3_BODY_POINTS.

if nargin < 3, dimmed = false; end

b = ch3_body_points(q);

if dimmed
    c_st = [0.85 0.6 0.6]; c_sw = [0.6 0.68 0.85]; c_k = [0.55 0.55 0.55];
    lw = 1.5;
else
    c_st = [0.85 0.2 0.2]; c_sw = [0.2 0.4 0.85]; c_k = [0 0 0];
    lw = 2.5;
end

plot(ax, [b.hip(1) b.torso_top(1)], [b.hip(2) b.torso_top(2)], '-', ...
     'Color', c_k, 'LineWidth', lw+0.5);
plot(ax, [b.hip(1) b.stance_knee(1) b.stance_foot(1)], ...
         [b.hip(2) b.stance_knee(2) b.stance_foot(2)], '-o', ...
     'Color', c_st, 'LineWidth', lw, 'MarkerSize', 5);
plot(ax, [b.hip(1) b.swing_knee(1) b.swing_foot(1)], ...
         [b.hip(2) b.swing_knee(2) b.swing_foot(2)], '--s', ...
     'Color', c_sw, 'LineWidth', lw-0.5, 'MarkerSize', 5);
plot(ax, b.hip(1), b.hip(2), 'o', 'Color', c_k, ...
     'MarkerFaceColor', c_k, 'MarkerSize', 6);
end
