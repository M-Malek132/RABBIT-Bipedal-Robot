function figs = ch5_plot_springmass(runs, savedir)
%CH5_PLOT_SPRINGMASS  Figure 5.3: the relative-degree-6 validation.
%
%   figs = ch5_plot_springmass(runs)
%   figs = ch5_plot_springmass(runs, savedir)
%
% runs is a struct array with .name, .sim, .label -- typically the CLF-QP
% baseline followed by the ECBF runs at two constraint levels, which is exactly
% the (a)/(b)/(c) of Fig. 5.3.
%
%   fig 1  Fig. 5.3 itself: cart positions, the constrained cart against its
%          limit, and the input force, one column per run
%   fig 2  the barrier record the thesis figure does not show -- h(t) and
%          y_rb(t) -- which is where the claim is actually checkable
%
% ------------------------------------------------------------- how to read fig 1
% The middle row is the whole figure. x3 (heavy) against x3d (dashed) and x3max
% (red). The baseline column crosses the red line; the ECBF columns do not, and
% the thesis's point is that they do not while THE BOTTOM ROW STAYS THE SAME
% SHAPE AND SCALE across columns. Fig. 5.3's caption makes that explicit --
% "varying the safety constraint while keeping the poles fixed keeps the peak
% forces and speed of system response the same" -- so the force axis is SHARED
% across columns. Letting each column autoscale would hide the one thing the
% figure is claiming.
%
% ---------------------------------------------------------------- and fig 2
% h(t) >= 0 is the claim; y_rb(t) >= 0 is the row the QP actually enforces
% (Remark 5.6). Plotting both separates "the constraint held" from "the
% constraint held BECAUSE the barrier was doing its job" -- on the baseline the
% first can be true by luck for a while, and the second is never even defined.
% y_rb is drawn on a symmetric log scale because it spans several decades as
% the row goes from slack to binding.
%
% See also CH5_PLOT_PENDULUM, CH5_REPORT, CH5_MAIN.

if nargin < 2, savedir = ''; end

n = numel(runs);
figs = gobjects(0);

%% ============================================================ fig 1: Fig 5.3
f1 = ch5_figure('ch5: Fig 5.3 serial spring-mass (relative degree 6)', 'page', 12.5);
t1 = tiledlayout(f1, 3, n, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(3, n);

u_lim = ch5_robust_ulim(runs);

x_lim = [inf -inf];
for i = 1:n
    x_lim(1) = min(x_lim(1), min(runs(i).sim.x(1:3,:), [], 'all'));
    x_lim(2) = max(x_lim(2), max(runs(i).sim.x(1:3,:), [], 'all'));
end
x_lim = x_lim + [-0.1 0.1]*diff(x_lim);

for i = 1:n
    s   = runs(i).sim;
    lvl = s.p.constraint.value;
    x3d = s.p.plant.x3d;

    % --- row 1: all three carts
    ax(1,i) = nexttile(t1, i); hold on; grid on;
    plot(s.t, s.x(1,:), 'LineWidth', 0.9);
    plot(s.t, s.x(2,:), 'LineWidth', 0.9);
    plot(s.t, s.x(3,:), 'LineWidth', 1.3);
    ylim(x_lim); xlim([0 s.t(end)]);
    title(run_title(runs(i).label));
    if i == 1
        ylabel('cart positions (m)');
        % two columns, in the corner the carts have left by t = 8 s
        legend({'x_1','x_2','x_3'}, 'Location','southeast', 'NumColumns', 2);
    end

    % --- row 2: the constrained cart against its limit
    %
    % xlim is set EXPLICITLY on every panel. Left to autoscale, a column whose
    % trace leaves the zoomed y-range early gets a shorter time axis than its
    % neighbours, and a figure whose entire purpose is comparing columns then
    % compares them on different axes.
    ax(2,i) = nexttile(t1, n+i); hold on; grid on;
    plot(s.t, s.x(3,:), 'LineWidth', 1.3);
    yline(x3d, 'k--', 'LineWidth', 0.8);
    yline(lvl, 'r-',  'LineWidth', 1.1);
    ylim([min(2.0, lvl-0.2), max([lvl, max(s.x(3,:))]) + 0.15]);
    xlim([0 s.t(end)]);
    if i == 1
        ylabel('x_3 (m)');
        % two columns, below the settling x_3 and right of its rise: stacked
        % in one, the entries reached up into the first; in one row, across
        % the second
        legend({'x_3','x_{3d}','x_3^{max}'}, 'Location','southeast', 'NumColumns', 2);
    end

    % --- row 3: input force, SHARED axis (see the header)
    ax(3,i) = nexttile(t1, 2*n+i); hold on; grid on;
    plot(s.t, s.u(1,:), 'LineWidth', 1.0, 'Color', [0.85 0.33 0.10]);
    ylim([-u_lim u_lim]);
    xlim([0 s.t(end)]);
    ch5_note_clipping(ax(3,i), s, u_lim, 'N');
    if i == 1, ylabel('u (N)'); end
end
xlabel(t1, 'Time (s)');
for i = 1:n, share_time(ax(:,i)); end    % a column is one run

figs(end+1) = f1;

%% ====================================================== fig 2: barrier record
f2 = ch5_figure('ch5: springmass barrier record', 'page', 10);
t2 = tiledlayout(f2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax2 = gobjects(2, 1);
cols = lines(max(n, 3));

ax2(1) = nexttile(t2); hold on; grid on;
for i = 1:n
    plot(runs(i).sim.t, runs(i).sim.h, 'Color', cols(i,:), 'LineWidth', 1.1, ...
         'DisplayName', runs(i).label);
end
yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylabel('h (m)');
title('The claim: h = x_3^{max} - x_3 \geq 0 for all t');

% Each run in ITS OWN colour. Left to the colour order, skipping the baseline
% (which has no barrier row) drew every ECBF run here in the colour the panel
% above gives the run before it.
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
    names = {'ch5_springmass_fig53', 'ch5_springmass_barrier'};
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(savedir, [names{k} '.png']), ...
                       'Resolution', 600);
    end
end

end

% ---------------------------------------------------------------------------
function c = run_title(label)
%RUN_TITLE  A run label as a column title, a line per comma-separated part.
% On one line "ECBF-CLF-QP, x_3^{max}-x_{3d} = 15 cm" is wider than a column
% at print size. The labels are TeX, which is the title's default interpreter.
c = strtrim(strsplit(label, ','));
end

function share_time(ax)
%SHARE_TIME  Stacked panels on one time axis, labelled on the bottom panel.
% The ticks are copied from the bottom panel rather than left to each axes: an
% axes whose tick labels are hidden picks its ticks differently, and its grid
% lines then disagree with the panel below. drawnow first, since the ticks
% MATLAB picks depend on the size the layout leaves the panel.
linkaxes(ax, 'x');
drawnow;
set(ax, 'XTick', ax(end).XTick);
set(ax(1:end-1), 'XTickLabel', {});
end
