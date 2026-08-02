function make_yd_smoothness_figure()
%MAKE_YD_SMOOTHNESS_FIGURE  Regenerate docs/figures/fig_yd_smoothness.png.
%
%   make_yd_smoothness_figure()
%
% Documents the smoothness of ch3_yd on both bases, for the verification
% section of RABBIT_documentation.tex.  Four panels:
%
%   (a) the four Bezier output profiles yd(s), with the nominal [0,1] band
%   (b) yd for one output, Bezier vs B-spline -- the B-spline CLAMPS s to
%       [0,1], so its value freezes outside the band while the Bezier
%       extrapolates smoothly
%   (c) dyd for the same output -- the B-spline value freezes but its
%       reported SLOPE freezes at the boundary slope rather than going to
%       zero, so outside [0,1] the returned dyd is not the derivative of the
%       returned yd
%   (d) |d2yd_bspline - d2yd_analytic| -- the B-spline second derivative is
%       central-differenced, and the difference is stepped in from the ends,
%       so accuracy collapses at s = 1, which is touchdown on every step
%
% Alpha comes from Results/ch3_reference_gait.mat.  Re-run after that gait
% changes, then recompile the document.
%
% See also CH3_YD, CH3_BEZIER, MAKE_DOC_FIGURES.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

figdir = fullfile(root, 'docs', 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end

S = load(fullfile(root, 'Results', 'ch3_reference_gait.mat'));
[~, ~, alpha] = ch3_col_unpack(S.z_opt, S.p);

pb = S.p;  pb.basis = 'bezier';
ps = S.p;  ps.basis = 'bspline';  ps.bsp_deg = S.p.bez_deg;   % equivalent curve

io   = 1;                                  % output shown in (b)-(d)
name = 'q_1 (stance hip)';

sw  = linspace(-0.35, 1.35, 900);          % wide, to show the clamp

% The central difference in the B-spline branch uses h = 1e-5 and clips its
% stencil to [0,1], so the one-sided degradation lives within one step of the
% endpoint.  Sample the DISTANCE from s = 1 logarithmically or it is invisible.
gap  = logspace(-8, -0.05, 700);
sin_ = 1 - gap;

n   = size(alpha,2) - 1;
deg = ps.bsp_deg;

YB = zeros(4, numel(sw));  DB = YB;
YS = zeros(1, numel(sw));  DS = YS;  DSold = YS;
for k = 1:numel(sw)
    [YB(:,k), DB(:,k)] = ch3_yd(alpha, sw(k), pb);
    [a, b]             = ch3_yd(alpha, sw(k), ps);
    YS(k) = a(io);  DS(k) = b(io);
    DSold(k) = former_dyd(alpha, n, deg, sw(k), io);
end

err = zeros(1, numel(sin_));  errOld = err;
for k = 1:numel(sin_)
    [~,~,a2] = ch3_yd(alpha, sin_(k), pb);
    [~,~,b2] = ch3_yd(alpha, sin_(k), ps);
    err(k)    = abs(b2(io) - a2(io));
    errOld(k) = abs(former_d2yd(alpha, n, deg, sin_(k), io) - a2(io));
end

cB = [0.15 0.35 0.75];      % bezier
cS = [0.85 0.33 0.10];      % bspline, current
cO = [0.72 0.72 0.72];      % bspline, former behaviour (for contrast)
band = [0.93 0.93 0.93];

f = figure('Color','w','Visible','off','Position',[100 100 1000 720]);

% ---- (a) the four profiles ------------------------------------------------
ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, sw, YB, 'LineWidth', 1.4);
xlabel(ax,'phase s'); ylabel(ax,'y_d(s)  [rad]');
title(ax,'(a) Bezier output profiles');
legend(ax, {'q_1','q_2','q_3','q_4'}, 'Location','northwest','Box','off');
xlim(ax,[sw(1) sw(end)]);  shade01(ax, band);

% ---- (b) value: clamp freeze ---------------------------------------------
ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, sw, YB(io,:), '-',  'Color', cB, 'LineWidth', 1.8);
plot(ax, sw, YS,       '--', 'Color', cS, 'LineWidth', 1.8);
xlabel(ax,'phase s'); ylabel(ax,'y_d(s)  [rad]');
title(ax, sprintf('(b) value, %s', name));
legend(ax, {'bezier','bspline (clamped)'}, 'Location','northwest','Box','off');
xlim(ax,[sw(1) sw(end)]);  shade01(ax, band);

% ---- (c) slope: consistent with the clamped value -------------------------
ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, sw, DSold,    '-',  'Color', cO, 'LineWidth', 1.4);
plot(ax, sw, DB(io,:), '-',  'Color', cB, 'LineWidth', 1.8);
plot(ax, sw, DS,       '--', 'Color', cS, 'LineWidth', 1.8);
yline(ax, 0, ':', 'Color',[0.4 0.4 0.4], 'HandleVisibility','off');
xlabel(ax,'phase s'); ylabel(ax,'dy_d/ds');
title(ax,'(c) slope: zero where the value is frozen');
legend(ax, {'bspline, former','bezier','bspline, current'}, ...
       'Location','northwest','Box','off');
xlim(ax,[sw(1) sw(end)]);  shade01(ax, band);

% ---- (d) second-derivative error near touchdown ---------------------------
ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
loglog(ax, gap, max(errOld,1e-16), '-', 'Color', cO, 'LineWidth', 1.4);
loglog(ax, gap, max(err,   1e-16), '-', 'Color', cS, 'LineWidth', 1.6);
set(ax,'XScale','log','YScale','log','XDir','reverse');
xline(ax, 1e-5, '--', 'Color',[0.4 0.4 0.4], 'HandleVisibility','off');
xlabel(ax,'distance from touchdown,  1 - s'); ylabel(ax,'|d^2y_d/ds^2 error|');
title(ax,'(d) bspline d^2y_d vs analytic');
legend(ax, {'clipped stencil (former)','interior stencil (current)'}, ...
       'Location','northwest','Box','off');
text(ax, 1.4e-5, 3e-11, 'FD step h = 10^{-5}', ...
     'FontSize',9, 'Color',[0.3 0.3 0.3], 'HorizontalAlignment','right');
xlim(ax,[min(gap) max(gap)]);

exportgraphics(f, fullfile(figdir,'fig_yd_smoothness.png'), 'Resolution', 150);
close(f);
fprintf('wrote %s\n', fullfile(figdir,'fig_yd_smoothness.png'));

end

% ---------------------------------------------------------------------------
function d = former_dyd(alpha, n, deg, s, io)
%FORMER_DYD  The pre-fix behaviour: evaluate the slope at the CLAMPED s, so it
% freezes at the boundary value instead of going to zero. Kept only to draw the
% contrast in panel (c).
sc = min(max(s, 0), 1);
v  = alpha * BSpline_derivative(n, deg, sc).';
d  = v(io);
end

% ---------------------------------------------------------------------------
function d = former_d2yd(alpha, n, deg, s, io)
%FORMER_D2YD  The pre-fix second derivative: a central difference whose stencil
% is CLIPPED against [0,1], so it silently degenerates to a half-width one-sided
% formula near the endpoints. Kept only to draw the contrast in panel (d).
h  = 1e-5;
sm = min(max(s-h, 0), 1);
sp = min(max(s+h, 0), 1);
if sp == sm
    d = 0;
else
    v = (alpha * BSpline_derivative(n, deg, sp).' - ...
         alpha * BSpline_derivative(n, deg, sm).') / (sp - sm);
    d = v(io);
end
end

% ---------------------------------------------------------------------------
function shade01(ax, c)
%SHADE01  Shade the nominal phase interval [0,1] behind the data.
%
% Call this AFTER the data is plotted and the limits are settled: the patch is
% sized to the current ylim and then pushed to the bottom of the stack. Drawing
% it first, with a nominally "large" height, makes the patch itself set the
% y-limits and flattens every curve in the axes.
yl = ylim(ax);
h  = patch(ax, [0 1 1 0], [yl(1) yl(1) yl(2) yl(2)], c, ...
           'EdgeColor','none', 'HandleVisibility','off');
uistack(h, 'bottom');
ylim(ax, yl);
end
