function figs = ch5_plot_pendulum(runs, savedir)
%CH5_PLOT_PENDULUM  Figure 5.4: the relative-degree-4 nonlinear validation.
%
%   figs = ch5_plot_pendulum(runs)
%   figs = ch5_plot_pendulum(runs, savedir)
%
%   fig 1  Fig. 5.4 itself, one column per run:
%            row 1  control outputs y1, y2
%            row 2  control inputs tau1, tau2
%            row 3  the end effector's path in the plane, against p2min
%   fig 2  the barrier record: py2(t) against its limit, and y_rb(t)
%
% ---------------------------------------------------------- how to read fig 1
% ROW 3 IS THE FIGURE. The workspace trace shows where the end effector was
% allowed to go, and the red line is where it was not. In the baseline column
% the trace sweeps a full circle of radius 2 and dives to py = -2; in the ECBF
% columns it flattens against the red line and slides along it. That sliding is
% the barrier being ACTIVE rather than merely satisfied.
%
% ROW 1 EXPLAINS THE COST. y2 = theta2 - theta2d is what the CLF wants at zero,
% and in the ECBF columns it is driven to roughly -2 rad in the middle of the
% run -- the arm folding. The CLF is not being satisfied there and it is not
% supposed to be: the barrier row is hard, the CLF row is slacked, and this is
% what that asymmetry looks like in a trajectory. Comparing the columns of row
% 1 without row 3 next to them would read as the controller doing worse, when
% it is doing something different and harder.
%
% The torque axis in row 2 is shared across columns, for the same reason as in
% ch5_plot_springmass: cross-column comparison is the point.
%
% See also CH5_PLOT_SPRINGMASS, CH5_ANIMATE_PENDULUM, CH5_REPORT.

if nargin < 2, savedir = ''; end

n = numel(runs);
figs = gobjects(0);

%% ============================================================ fig 1: Fig 5.4
% LAID OUT BY HAND, in centimetres, not with tiledlayout. tiledlayout gives
% every row room for the tallest title and the tallest axis label of ANY row
% -- here the two-line column titles and row 2's time label -- and with the
% square workspace panels that left the two time rows 1.3 cm tall between 2 cm
% gaps. Each margin below holds only what sits in it.
L = struct('left', 1.6, 'right', 0.2, 'gapx', 0.85, ...  % y labels | tick labels
           'top', 1.05, ...          % the two-line column titles
           'row', 2.6, ...           % height of each time row
           'note', 0.85, ...         % row 2's two-line clipping note
           'time', 1.15, ...         % row 2's ticks and time label
           'bottom', 1.75);          % row 3's ticks and label, then the legend
w  = (17.4 - L.left - L.right - (n-1)*L.gapx) / n;   % column; row 3 is w x w
y3 = L.bottom;
y2 = y3 + w + L.time;
y1 = y2 + L.row + L.note;

f1 = ch5_figure('ch5: Fig 5.4 two-link pendulum, elastic actuators', 'page', ...
                y1 + L.row + L.top);
ax = gobjects(3, n);

u_lim = ch5_robust_ulim(runs);

for i = 1:n
    s    = runs(i).sim;
    thd  = s.p.plant.thetad(:);
    lvl  = s.p.constraint.value;
    reach = sum(s.p.plant.l);
    x0   = L.left + (i-1)*(w + L.gapx);

    % --- row 1: control outputs
    ax(1,i) = axes(f1, 'Units', 'centimeters', 'Position', [x0 y1 w L.row]);
    hold on; grid on;
    plot(s.t, s.x(1,:) - thd(1), 'LineWidth', 1.0);
    plot(s.t, s.x(2,:) - thd(2), 'LineWidth', 1.0);
    yline(0, 'k:');
    xlim([0 s.t(end)]);
    title(run_title(runs(i).label));
    if i == 1
        ylabel('outputs y (rad)');     % short enough for a 2.6 cm row
        legend({'y_1 = \theta_1 - \theta_{1d}', 'y_2 = \theta_2 - \theta_{2d}'}, ...
               'Location', 'southeast');
    end

    % --- row 2: control inputs, SHARED axis
    %
    % Drawn thin because the ECBF traces CHATTER. That is not noise to be
    % filtered out of the plot: while the barrier row is active the controller
    % is holding the state on the constraint boundary with a zero-order hold,
    % the state drifts off between samples, and the next sample corrects it --
    % a discrete-time sliding mode. It is the honest picture of what enforcing
    % a hard constraint at a finite rate looks like, and the %act column in
    % ch5_report says how much of the run is spent doing it.
    ax(2,i) = axes(f1, 'Units', 'centimeters', 'Position', [x0 y2 w L.row]);
    hold on; grid on;
    plot(s.t, s.u(1,:), 'LineWidth', 0.5);
    plot(s.t, s.u(2,:), 'LineWidth', 0.5);
    ylim([-u_lim u_lim]);
    xlim([0 s.t(end)]);
    ch5_note_clipping(ax(2,i), s, u_lim, 'Nm');
    xlabel('Time (s)');
    if i == 1
        ylabel('torques \tau (Nm)');
        % top right: clear of the first column's torques after t = 11 s
        legend({'\tau_1','\tau_2'}, 'Location','northeast');
    end

    % --- row 3: the end effector in the plane
    ax(3,i) = axes(f1, 'Units', 'centimeters', 'Position', [x0 y3 w w]);
    hold on; grid on; axis equal;
    P   = ch5_pend_points(s.x(1:2,:), s.p);
    tip = squeeze(P(:,3,:));

    % the reachable circle, for scale
    a = linspace(0, 2*pi, 200);
    plot(reach*cos(a), reach*sin(a), ':', 'Color', [0.6 0.6 0.6]);

    plot(tip(1,:), tip(2,:), 'LineWidth', 1.0, 'Color', [0 0.45 0.74]);
    yline(lvl, 'r-', 'LineWidth', 1.1);
    plot(tip(1,1),   tip(2,1),   'ko', 'MarkerFaceColor','g', 'MarkerSize', 5);
    plot(tip(1,end), tip(2,end), 'ks', 'MarkerFaceColor','r', 'MarkerSize', 5);

    % the arm at the instant the constraint was tightest -- the pose the
    % barrier forced, which is otherwise invisible in a trace of the tip alone
    [~, k] = min(s.h);
    plot(squeeze(P(1,:,k)), squeeze(P(2,:,k)), '-o', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'MarkerSize', 3, ...
         'MarkerFaceColor', [0.2 0.2 0.2]);

    xlim([-1.1 1.1]*reach); ylim([-1.1 1.1]*reach);
    xlabel('x (m)');
    if i == 1, ylabel('y (m)'); end
end
% Centred below the figure rather than under the first panel, where it shrank
% that panel alone and the three workspaces came out at different sizes.
lg = legend(ax(3,1), {'reach','p_2 path','p_{2min}','start','end','pose at min h'}, ...
            'NumColumns', 6);
lg.Units = 'centimeters';
lg.Position(1:2) = [(17.4 - lg.Position(3))/2, 0.05];
for i = 1:n, share_time(ax(1:2,i)); end  % the workspace row has no time axis

figs(end+1) = f1;

%% ====================================================== fig 2: barrier record
f2 = ch5_figure('ch5: pendulum barrier record', 'page', 10);
t2 = tiledlayout(f2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax2 = gobjects(2, 1);
cols = lines(max(n, 3));

ax2(1) = nexttile(t2); hold on; grid on;
for i = 1:n
    s = runs(i).sim;
    plot(s.t, s.h + s.p.constraint.value, 'Color', cols(i,:), 'LineWidth', 1.1, ...
         'DisplayName', runs(i).label);
end
% each limit once: the baseline runs against the first ECBF run's level
for lvl = unique(arrayfun(@(r) r.sim.p.constraint.value, runs))
    yline(lvl, 'r--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
ylabel('p^y_2 (m)');
title('End-effector height against its limits');
% a tick every metre, so the limits can be read off; left alone the axis
% labelled only -2, 0 and 2
yl = ylim(ax2(1));
yticks(ax2(1), ceil(yl(1)):floor(yl(2)));

% each run in its own colour, as in ch5_plot_springmass: left to the colour
% order, skipping the baseline shifted every ECBF run to the wrong colour
ax2(2) = nexttile(t2); hold on; grid on;
for i = 1:n
    y = runs(i).sim.y_rb;
    if all(isnan(y)), continue; end
    plot(runs(i).sim.t, sign(y).*log10(1 + abs(y)), 'Color', cols(i,:), ...
         'LineWidth', 1.1);
end
yline(0, 'k-', 'LineWidth', 0.8);
ylabel({'$\mathrm{sgn}(y_{r_b})\,\times$', '$\log_{10}(1+|y_{r_b}|)$'}, ...
       'Interpreter', 'latex');
title('The row the QP enforces (Remark 5.6): y_{r_b} \geq 0; zero means binding');
xlabel(t2, 'Time (s)');
% one legend for both panels, above them; TeX, since the labels are written in it
lg = legend(ax2(1), 'Orientation', 'horizontal');
lg.Layout.Tile = 'north';
share_time(ax2);

figs(end+1) = f2;

%% ------------------------------------------------------------------- save
if ~isempty(savedir)
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    names = {'ch5_pendulum_fig54', 'ch5_pendulum_barrier'};
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(savedir, [names{k} '.png']), ...
                       'Resolution', 600);
    end
end

end

% ---------------------------------------------------------------------------
function c = run_title(label)
%RUN_TITLE  A run label as a column title, a line per comma-separated part.
% On one line "ECBF-CLF-QP, p_{2min} = -1.0 m" crowds a column at print size.
% The labels are TeX, which is the title's default interpreter.
c = strtrim(strsplit(label, ','));
end

function share_time(ax)
%SHARE_TIME  Stacked panels on one time axis, labelled on the bottom panel.
% As in ch5_plot_springmass: ticks copied from the bottom panel, after drawnow,
% so a panel whose labels are hidden keeps the same grid.
linkaxes(ax, 'x');
drawnow;
set(ax, 'XTick', ax(end).XTick);
set(ax(1:end-1), 'XTickLabel', {});
end
