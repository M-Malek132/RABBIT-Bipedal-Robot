function ch3_doc_figures()
%CH3_DOC_FIGURES  Redraw the Chapter-3 report figures into docs/figures.
%
%   ch3_doc_figures()
%
% The Persian report (docs/ch3_report_fa.tex) shows three figures. The pipeline
% diagram is TikZ inside the .tex; the other two are drawn here. Then rebuild:
%
%   ch3_doc_figures
%   !cd docs && latexmk -xelatex ch3_report_fa.tex
%
%   docs/figures/              drawn from                              .tex width
%   ch3_gait_posture195.png    Results/ch3_gait_posture_195.mat         \textwidth       17.4 cm
%   ch3_defect_order.png       seed rollouts, six meshes, two bases     0.82\textwidth   14.3 cm
%
% DRAWN AT PRINT SIZE, TEXT IN POINTS. LaTeX scales a figure to the width it
% sets, so whether its text can be read depends on the text's size relative to
% the figure, not on the resolution. Exported from ch3_plot_gait's 1250-pixel
% window at MATLAB's default font size, the six gait panels printed at 3-4 pt
% across the page. Drawn at 17.4 cm with 9 pt text, they print at 9 pt.
%
% THE DEFECT STUDY IS RECOMPUTED, NOT LOADED. No result file holds it, and its
% twelve seed rollouts take seconds. It is the report's order study: the
% degree-3 B-spline default against the Bezier basis (the degree-5 B-spline is
% the same curve as the Bezier and measures identically), at N = 9 ... 65. The
% numbers are printed, so they can be checked against the report's table.
%
% Requires a normal (JVM-enabled) MATLAB session.
%
% See also CH3_PLOT_GAIT, CH3_COL_SEED, CH3_COL_EVAL.

root   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(root, 'docs', 'figures');

%% the reference gait, six panels across the page
S = load(fullfile(root, 'Results', 'ch3_gait_posture_195.mat'), 'z', 'p');
p = ch3_upgrade_params(S.p);
E = ch3_col_eval(S.z, p);
[X, T, alpha] = ch3_col_unpack(S.z, p);
N = E.N;
t = linspace(0, T, N);

f  = print_figure(17.4, 11.4);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. stick figure
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
idx  = round(linspace(1, N, 6));
cols = parula(numel(idx) + 1);
for j = 1:numel(idx)
    b = ch3_body_points(X(1:p.nq, idx(j)));
    c = cols(j,:);
    plot(ax, [b.hip(1) b.torso_top(1)], [b.hip(2) b.torso_top(2)], '-', 'Color', c, 'LineWidth', 1.4);
    plot(ax, [b.hip(1) b.stance_knee(1) b.stance_foot(1)], ...
             [b.hip(2) b.stance_knee(2) b.stance_foot(2)], '-o', 'Color', c, 'MarkerSize', 2);
    plot(ax, [b.hip(1) b.swing_knee(1) b.swing_foot(1)], ...
             [b.hip(2) b.swing_knee(2) b.swing_foot(2)], '--s', 'Color', c, 'MarkerSize', 2);
end
yline(ax, 0, 'k-', 'LineWidth', 0.8);
title(ax, sprintf('%.3f m in %.3f s', E.L_step, T));
xlabel(ax, 'x [m]'); ylabel(ax, 'z [m]');

% 2. virtual constraints
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
s_nodes = zeros(1, N); Y = zeros(p.ny, N); YD = zeros(p.ny, N);
for k = 1:N
    s_nodes(k) = ch3_phase(X(:,k), p);
    Y(:,k)  = p.H * X(1:p.nq, k);
    YD(:,k) = ch3_yd(alpha, s_nodes(k), p);
end
co = lines(p.ny);
for i = 1:p.ny
    plot(ax, s_nodes, Y(i,:),  '-',  'Color', co(i,:), 'LineWidth', 1.0);
    plot(ax, s_nodes, YD(i,:), '--', 'Color', co(i,:), 'LineWidth', 0.8);
end
xlabel(ax, 'phase s'); ylabel(ax, 'angle [rad]');
title(ax, 'q (solid), y_d(s) (dashed)');

% 3. torques
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
plot(ax, t, E.u.', 'LineWidth', 0.9);
if p.limits.enable.torque
    yline(ax,  p.limits.u_max, 'r--'); yline(ax, -p.limits.u_max, 'r--');
end
xlabel(ax, 't [s]'); ylabel(ax, 'u [Nm]');
title(ax, sprintf('torques, peak %.1f Nm', max(abs(E.u(:)))));
% One joint legend serves panels 2 and 3, outside both: output i is q(3+i), and
% u_i drives the same joint. Inside a 5 cm panel a legend covers the curves.
lgj = legend(ax, {'stance hip', 'stance knee', 'swing hip', 'swing knee'}, 'NumColumns', 4);
lgj.Layout.Tile = 'north';

% 4. ground reaction force
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
plot(ax, t, E.lam(1,:), 'LineWidth', 0.9);
plot(ax, t, E.lam(2,:), 'LineWidth', 0.9);
plot(ax, t,  p.limits.mu_s * E.lam(2,:), 'k--', 'LineWidth', 0.6);
plot(ax, t, -p.limits.mu_s * E.lam(2,:), 'k--', 'LineWidth', 0.6, 'HandleVisibility', 'off');
yline(ax, p.limits.Fz_min, 'r:', 'HandleVisibility', 'off');
yline(ax, 0, 'k-', 'HandleVisibility', 'off');
xlabel(ax, 't [s]'); ylabel(ax, 'force [N]');
title(ax, sprintf('GRF, min F_z %.0f N', min(E.lam(2,:))));
lgf = legend(ax, {'ground force F_x', 'F_z', 'friction cone \pm\mu_s F_z'}, 'NumColumns', 3);
lgf.Layout.Tile = 'south';

% 5. distance from Z
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
eta = zeros(2*p.ny, N);
for k = 1:N
    [yk, ydk] = ch3_outputs(X(:,k), alpha, p);
    eta(:,k) = [yk; ydk];
end
plot(ax, t, max(abs(eta), [], 1), 'LineWidth', 1.0);
set(ax, 'YScale', 'log', 'YLim', [1e-8 1e-4], 'YTick', 10.^(-8:2:-4));
xlabel(ax, 't [s]'); ylabel(ax, 'max |\eta|');
title(ax, 'distance from Z');

% 6. swing-foot height
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
plot(ax, t, E.sw_h, 'LineWidth', 1.0);
yline(ax, p.limits.clearance, 'r--');
yline(ax, 0, 'k-');
xlabel(ax, 't [s]'); ylabel(ax, 'height [m]');
title(ax, 'swing-foot height');

save_as(f, figdir, 'ch3_gait_posture195.png');
close(f);

%% Hermite-Simpson defect under mesh refinement, on the seed rollout
Ns    = [9 13 21 33 49 65];
bases = {'bspline', 3; 'bezier', 5};
[h, dmax, dmed] = deal(zeros(size(bases, 1), numel(Ns)));
for b = 1:size(bases, 1)
    for k = 1:numel(Ns)
        pn = ch3_params('N_nodes', Ns(k), 'basis', bases{b,1}, 'bsp_deg', bases{b,2});
        En = ch3_col_eval(ch3_col_seed(pn), pn);
        D  = max(abs(En.defect), [], 1);            % worst state row per interval
        h(b,k)    = En.h;
        dmax(b,k) = max(D);
        dmed(b,k) = median(D);
    end
end
order = log(dmax(:, 2:end) ./ dmax(:, 1:end-1)) ./ log(h(:, 2:end) ./ h(:, 1:end-1));
fprintf(' ch3_doc_figures: defect study (max | median)\n');
for k = 1:numel(Ns)
    fprintf('   N = %2d   B-spline 3: h %.4f  %.3e | %.3e    Bezier 5: h %.4f  %.3e | %.3e\n', ...
            Ns(k), h(1,k), dmax(1,k), dmed(1,k), h(2,k), dmax(2,k), dmed(2,k));
end

red  = [0.70 0.16 0.12];
blue = [0.13 0.36 0.65];
f  = print_figure(14.3, 7.8);
ax = axes(f); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
plot(ax, h(1,:), dmax(1,:), '-o',  'Color', red,  'MarkerFaceColor', red,  'MarkerSize', 3.5, 'LineWidth', 1.0);
plot(ax, h(1,:), dmed(1,:), '--o', 'Color', red,  'MarkerFaceColor', 'w',  'MarkerSize', 3.5, 'LineWidth', 0.9);
plot(ax, h(2,:), dmax(2,:), '-o',  'Color', blue, 'MarkerFaceColor', blue, 'MarkerSize', 3.5, 'LineWidth', 1.0);
plot(ax, h(2,:), dmed(2,:), '--o', 'Color', blue, 'MarkerFaceColor', 'w',  'MarkerSize', 3.5, 'LineWidth', 0.9);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.0095 0.095], 'YLim', [1e-9 1e-1], ...
        'XTick', [0.01 0.02 0.05], 'YTick', 10.^(-9:2:-1));
ax.XAxis.MinorTickValues = [0.03 0.04 0.06 0.07 0.08 0.09];
xlabel(ax, sprintf('mesh interval h [s]   (N = %d on the left, N = %d on the right)', Ns(end), Ns(1)));
ylabel(ax, 'Hermite–Simpson defect');
% The stall is where the B-spline's interior knots land inside an interval:
% the maximum's order between N = 13, 21 and 33.
text(ax, 0.0105, 3e-2, sprintf('maximum stalls: order %.2f, %.2f', order(1,2), order(1,3)), ...
     'Color', red, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
legend(ax, {'B-spline degree 3 (default), maximum', 'B-spline degree 3, median over intervals', ...
            'Bézier degree 5, maximum', 'Bézier degree 5, median over intervals'}, ...
       'Location', 'northoutside', 'NumColumns', 2);
save_as(f, figdir, 'ch3_defect_order.png');
close(f);
end

% ---------------------------------------------------------------------------
function f = print_figure(width_cm, height_cm)
%PRINT_FIGURE  A figure the size the report prints it, with 9 pt text. Legends
% carry no box, whose white fill hides the curves under it.
f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 width_cm height_cm], ...
           'DefaultAxesFontSize', 9, ...
           'DefaultAxesLabelFontSizeMultiplier', 1.1, ...
           'DefaultAxesTitleFontSizeMultiplier', 1.1, ...
           'DefaultAxesTitleFontWeight', 'normal', ...
           'DefaultTextFontSize', 9, ...
           'DefaultAxesLineWidth', 0.5, ...
           'DefaultLineLineWidth', 0.8, ...
           'DefaultLegendBox', 'off', ...
           'DefaultLegendItemTokenSize', [15 8]);
end

function save_as(f, figdir, name)
% exportgraphics crops to the content; saveas would keep the window's margins.
exportgraphics(f, fullfile(figdir, name), 'Resolution', 600);
fprintf(' ch3_doc_figures: %s\n', name);
end
