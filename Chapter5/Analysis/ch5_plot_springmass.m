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
f1 = figure('Name', 'ch5: Fig 5.3 serial spring-mass (relative degree 6)', ...
            'Position', [60 60 380*n 720]);

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
    subplot(3, n, i); hold on; grid on;
    plot(s.t, s.x(1,:), 'LineWidth', 1.2);
    plot(s.t, s.x(2,:), 'LineWidth', 1.2);
    plot(s.t, s.x(3,:), 'LineWidth', 1.8);
    ylim(x_lim); xlim([0 s.t(end)]);
    ylabel('cart positions  [m]');
    title(runs(i).label, 'Interpreter', 'none', 'FontWeight', 'normal');
    if i == 1, legend({'x_1','x_2','x_3'}, 'Location','southeast'); end

    % --- row 2: the constrained cart against its limit
    %
    % xlim is set EXPLICITLY on every panel. Left to autoscale, a column whose
    % trace leaves the zoomed y-range early gets a shorter time axis than its
    % neighbours, and a figure whose entire purpose is comparing columns then
    % compares them on different axes.
    subplot(3, n, n+i); hold on; grid on;
    plot(s.t, s.x(3,:), 'LineWidth', 1.8);
    yline(x3d, 'k--', 'LineWidth', 1.0);
    yline(lvl, 'r-',  'LineWidth', 1.5);
    ylim([min(2.0, lvl-0.2), max([lvl, max(s.x(3,:))]) + 0.15]);
    xlim([0 s.t(end)]);
    ylabel('x_3  [m]');
    if i == 1, legend({'x_3','x_{3d}','x_3^{max}'}, 'Location','southeast'); end

    % --- row 3: input force, SHARED axis (see the header)
    ax3 = subplot(3, n, 2*n+i); hold on; grid on;
    plot(s.t, s.u(1,:), 'LineWidth', 1.4, 'Color', [0.85 0.33 0.10]);
    ylim([-u_lim u_lim]);
    xlim([0 s.t(end)]);
    ch5_note_clipping(ax3, s, u_lim, 'N');
    xlabel('Time (s)'); ylabel('u  [N]');
end

figs(end+1) = f1;

%% ====================================================== fig 2: barrier record
f2 = figure('Name', 'ch5: springmass barrier record', 'Position', [80 80 900 520]);

subplot(2,1,1); hold on; grid on;
for i = 1:n
    plot(runs(i).sim.t, runs(i).sim.h, 'LineWidth', 1.5);
end
yline(0, 'k-', 'LineWidth', 1.2);
ylabel('h(x) = x_3^{max} - x_3   [m]');
title('The claim: h \geq 0 for all t');
legend({runs.label}, 'Interpreter','none', 'Location','best');

subplot(2,1,2); hold on; grid on;
for i = 1:n
    y = runs(i).sim.y_rb;
    if all(isnan(y)), continue; end
    plot(runs(i).sim.t, sign(y).*log10(1 + abs(y)), 'LineWidth', 1.5);
end
yline(0, 'k-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('sgn(y_{r_b})\cdot log_{10}(1+|y_{r_b}|)');
title('The row the QP enforces (Remark 5.6): y_{r_b} \geq 0; zero means binding');

figs(end+1) = f2;

%% ------------------------------------------------------------------- save
if ~isempty(savedir)
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    names = {'ch5_springmass_fig53', 'ch5_springmass_barrier'};
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(savedir, [names{k} '.png']), ...
                       'Resolution', 150);
    end
end

end
