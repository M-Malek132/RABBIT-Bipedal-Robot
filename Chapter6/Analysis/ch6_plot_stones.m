function fig = ch6_plot_stones(sim, file, ttl)
%CH6_PLOT_STONES  Step length against the foothold window.  Figs. 6.5, 6.25(a).
%
%   fig = ch6_plot_stones(sim)
%   fig = ch6_plot_stones(sim, file, ttl)
%
% Two panels, because the chapter uses two different views of the same fact:
%
%   TOP   -- per step (Fig. 6.5). Blue bars are the foothold window
%            [l_min, l_max] "indicating the size of the footholds", red markers
%            the achieved l_s. A marker inside its bar is a placement that
%            landed on the stone.
%
%   BOTTOM -- against time (Figs. 6.10, 6.25a). The swing foot's horizontal
%            position l_f(t) through every step, with l_min(t) and l_max(t)
%            drawn as they move. This is the panel that makes Section 6.2
%            legible: on a moving stone the window the step ENDS at is not the
%            one it started with, and only a time axis shows that.
%
% The bottom panel restarts l_f at each impact because l_f is measured from the
% CURRENT stance foot -- the discontinuity at every step boundary is the leg
% swap, not a jump in the robot.
%
% Inputs
%   sim  : output of ch6_simulate
%   file : optional path to save a PNG
%   ttl  : optional title string
%
% Output
%   fig : figure handle
%
% See also CH6_SIMULATE, CH6_PLOT_BARRIERS, CH6_REPORT.

if nargin < 2, file = ''; end
if nargin < 3 || isempty(ttl), ttl = 'Dynamic walking on stepping stones'; end

n = sim.n_ok;
if n == 0
    error('ch6_plot_stones:empty', 'The run completed no steps: %s', sim.reason);
end

fig = figure('Color', 'w', 'Position', [80 80 980 700]);

%% -------------------------------------------------------- per step (Fig 6.5)
subplot(2,1,1); hold on; box on;

for k = 1:n
    % The window as a bar, drawn the way Fig. 6.5 draws it.
    plot([k k], [sim.l_min(k) sim.l_max(k)], '-', ...
         'Color', [0.20 0.40 0.85], 'LineWidth', 7);
end
h_hit  = plot(find( sim.hit), sim.l_s( sim.hit), 'o', 'MarkerSize', 7, ...
              'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', 'k');
h_miss = plot(find(~sim.hit), sim.l_s(~sim.hit), 'x', 'MarkerSize', 11, ...
              'Color', [0.85 0.15 0.15], 'LineWidth', 2.5);

xlabel('Step Number'); ylabel('Step Length [m]');
title(sprintf('%s  --  %d/%d placements on the stone', ttl, sum(sim.hit), n));
xlim([0.4, n+0.6]);
lg = {'foothold [l_{min}, l_{max}]', 'l_s on the stone', 'l_s missed'};
keep = [true, ~isempty(h_hit), ~isempty(h_miss)];
hs = [plot(nan, nan, '-', 'Color', [0.20 0.40 0.85], 'LineWidth', 7), ...
      h_hit, h_miss];
legend(hs(keep(1:numel(hs))), lg(keep(1:numel(hs))), 'Location', 'best');
grid on;

%% ------------------------------------------------ against time (Figs 6.10, 6.25a)
subplot(2,1,2); hold on; box on;

t_off = 0;
for k = 1:n
    lg_k = sim.steps(k).log;
    tt   = t_off + lg_k.t;
    plot(tt, lg_k.l_f,   '-',  'Color', [0.10 0.25 0.70], 'LineWidth', 2);
    plot(tt, lg_k.l_min, '--', 'Color', [0.90 0.45 0.10], 'LineWidth', 1.2);
    plot(tt, lg_k.l_max, '--', 'Color', [0.90 0.45 0.10], 'LineWidth', 1.2);
    plot(t_off + sim.steps(k).T, sim.l_s(k), 'o', 'MarkerSize', 6, ...
         'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', 'k');
    t_off = t_off + sim.steps(k).T;
end
xlabel('Time (s)'); ylabel('[m]');
legend({'l_f (swing foot)', 'l_{min}(t), l_{max}(t)'}, 'Location', 'best');
grid on;

if ~isempty(file)
    exportgraphics(fig, file, 'Resolution', 150);
    fprintf('[ch6_plot_stones] wrote %s\n', file);
end

end
