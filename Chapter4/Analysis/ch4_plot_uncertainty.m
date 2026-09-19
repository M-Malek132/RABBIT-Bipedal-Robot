function figs = ch4_plot_uncertainty(C, p, savedir, case_names)
%CH4_PLOT_UNCERTAINTY  The Chapter-4 figures, from a comparison sweep.
%
%   figs = ch4_plot_uncertainty(C, p)
%   figs = ch4_plot_uncertainty(C, p, savedir)
%   figs = ch4_plot_uncertainty(C, p, savedir, case_names)
%
% Takes the struct array from ch4_compare_controllers and draws the figures the
% chapter uses to make its case:
%
%   fig 1  CLF vs time, one panel per model perturbation      (Figs 4.2, 4.8)
%   fig 2  output tracking errors, cases across columns       (Fig 4.3)
%   fig 3  joint torques, cases across columns                (Figs 4.4, 4.9)
%   fig 4  torso phase portrait across cases                  (Figs 4.6, 4.10)
%   fig 5  the L1 estimator: theta_hat against the true theta (only when the
%          sweep contains an adaptive controller)
%
% On a Case IV sweep (preset 'case4'), figs 2 and 3 are Fig 4.5.
%
% HOW TO READ FIGURE 1, since it is the one the chapter leans on. The claim is
% NOT that the robust/adaptive curve is lowest -- on Case I it need not be, and
% Remark 4.7 says as much. The claim is that ONE CONTROLLER'S CURVE HAS THE
% SAME SHAPE IN ALL THREE PANELS while the baseline's degrades as the panels
% go down. So compare a single colour across panels, not the colours within a
% panel. The y-axes are therefore shared across panels on purpose; rescaling
% each panel to its own data would destroy exactly the comparison being made.
%
% Inputs
%   C       : struct array from ch4_compare_controllers (needs .traj, i.e.
%             opts.store_traj left on)
%   p          : parameter struct, for labels
%   savedir    : optional directory; figures are written as ch4_fig*.png
%   case_names : optional cell of the chapter's case numerals, one per scale in
%                sweep order. The default numbers the scales I, II, III by
%                position, which matches both default presets; a Case IV sweep
%                passes {'IV'}.
%
% Output
%   figs : vector of figure handles
%
% See also CH4_COMPARE_CONTROLLERS, CH4_REPORT, CH3_PLOT_GAIT.

if nargin < 3, savedir = ''; end

if isempty(C) || all(cellfun(@isempty, {C.traj}))
    error('ch4_plot_uncertainty:noTraj', ...
          ['The sweep carries no trajectories. Re-run ' ...
           'ch4_compare_controllers with opts.store_traj = true.']);
end

names  = unique({C.name},  'stable');
scales = unique([C.mass_scale], 'stable');
nC = numel(names); nS = numel(scales);

if nargin < 4 || isempty(case_names)
    case_names = arrayfun(@case_of_scale, scales, 'UniformOutput', false);
elseif numel(case_names) ~= nS
    error('ch4_plot_uncertainty:caseNames', ...
          'case_names has %d entries for a sweep over %d scales.', ...
          numel(case_names), nS);
end

cols = lines(max(nC,3));
figs = gobjects(0);

% the time panels share one axis, out to the longest run, so a run that
% ends early shows as one
ran  = C(~cellfun(@isempty, {C.traj}));
tmax = max(arrayfun(@(e) e.traj.t(end), ran));

% ---------------------------------------------------------------- fig 1: CLF
f1 = ch4_figure('ch4: CLF under model perturbation', 'column', 1.9 + 2.75*nS);
t1 = tiledlayout(f1, nS, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
ax = gobjects(1,nS);
% LOG AXIS, deliberately. A failing controller reaches V ~ 2e5 while a working
% one sits near 1, so on a shared linear axis -- and the axis MUST be shared, or
% the cross-case comparison this figure exists for is destroyed -- every
% informative curve collapses onto zero and the plot shows a single spike. Five
% decades of range is the finding, not an inconvenience, so use an axis that can
% show five decades.
%
% The floor is set from the data rather than left at realmin. Each run STARTS ON
% THE ORBIT, where V is zero to numerical precision, so an unbounded log axis
% devotes 300 decades to the first sample and squashes everything that matters
% into the top inch.
Vmax = 0;
for is = 1:nS
    for ic = 1:nC
        e = pick(C, names{ic}, scales(is));
        if isempty(e) || isempty(e.traj), continue; end
        Vv = e.traj.V;
        if all(isnan(Vv)), Vv = e.traj.eta_n.^2; end
        Vmax = max(Vmax, max(Vv));
    end
end
if ~isfinite(Vmax) || Vmax <= 0, Vmax = 1; end
Vfloor = Vmax * 1e-5;

for is = 1:nS
    ax(is) = nexttile(t1); hold on; grid on;
    for ic = 1:nC
        e = pick(C, names{ic}, scales(is));
        if isempty(e) || isempty(e.traj), continue; end
        Vv = e.traj.V;
        if all(isnan(Vv)), Vv = e.traj.eta_n.^2; end
        plot(e.traj.t, max(Vv, Vfloor), 'Color', cols(ic,:), ...
             'DisplayName', names{ic});
    end
    % no minor grid: at print size nine lines a decade bury the curves
    set(gca, 'YScale', 'log', 'YMinorGrid', 'off', 'YMinorTick', 'off');
    ylim([Vfloor, Vmax*2]);
    decade_ticks(gca);
    title(sprintf('Case %s: model scale = %.2g', case_names{is}, scales(is)));
end
ylabel(t1, '$V_\varepsilon$', 'Interpreter', 'latex');
xlabel(t1, 'Time (s)');
% Two short rows at the top of the first panel, which the sweeps leave empty:
% in a layout tile of its own a legend took a centimetre from the panels, one
% row across the panel covered Case IV's baselines, and one column reached down
% to Case I's footstrike spikes.
legend(ax(1), 'Interpreter', 'none', 'NumColumns', 2, 'Location', 'north');
share_time(ax, [0 tmax]);
figs(end+1) = f1;

% ------------------------------------------------- fig 2 / 3: outputs, torques
f2 = panel_grid(C, names, scales, cols, 'y',  'y_%d (deg)', 180/pi, ...
                'ch4: tracking errors');
f3 = panel_grid(C, names, scales, cols, 'u',  'u_%d (Nm)',  1, ...
                'ch4: joint torques');
figs(end+1) = f2;
figs(end+1) = f3;

% ------------------------------------------------------ fig 4: phase portrait
f4 = ch4_figure('ch4: torso phase portrait', 'column', 8.8);
t4 = tiledlayout(f4, 1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
ax4 = nexttile(t4); hold on; grid on;
mk = {'-','--',':','-.'};
for ic = 1:nC
    for is = 1:nS
        e = pick(C, names{ic}, scales(is));
        if isempty(e) || isempty(e.traj), continue; end
        qt  = e.traj.x(3,  :) * 180/pi;
        dqt = e.traj.x(3+p.nq, :) * 180/pi;
        plot(qt, dqt, mk{min(is,numel(mk))}, 'Color', cols(ic,:), ...
             'LineWidth', 0.6, 'HandleVisibility', 'off');
    end
end
% The legend keys colour to controller and line style to case instead of
% listing every run: nC x nS entries covered half the portrait.
for ic = 1:nC
    plot(NaN, NaN, '-', 'Color', cols(ic,:), 'LineWidth', 1.2, ...
         'DisplayName', names{ic});
end
if nS > 1
    for is = 1:nS
        plot(NaN, NaN, mk{min(is,numel(mk))}, 'Color', [0.25 0.25 0.25], ...
             'LineWidth', 1.2, 'DisplayName', sprintf('scale %.2g', scales(is)));
    end
end
xlabel('q_{torso} (deg)'); ylabel('dq_{torso} (deg/s)');
% NOT "uncertainty moves the periodic orbit", which is how Fig. 4.6 reads. A
% uniform mass scale leaves the hybrid zero dynamics exactly invariant (both
% halves are asserted in ch4_test_model), so a controller that tracks walks
% the same orbit in every case, and whatever separates the cases here is
% tracking error. Only a non-uniform change such as load_mass moves the orbit.
% (The report's caption says so; a subtitle saying it does not fit a column.)
title('Torso phase portrait');
legend(ax4, 'Interpreter', 'none', 'Orientation', 'horizontal', ...
       'NumColumns', max(nC, nS), 'Location', 'southoutside');
figs(end+1) = f4;

% ------------------------------------------------------- fig 5: L1 estimator
has_l1 = any(cellfun(@(n) any(strcmpi(n,{'l1','l1_con'})), names));
if has_l1
    f5 = ch4_figure('ch4: L1 estimator', 'column', 2.2 + 2.5*nS);
    t5 = tiledlayout(f5, nS, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    ax5 = gobjects(1, nS);
    est = false(1, nC);                       % controllers that estimate
    for is = 1:nS
        ax5(is) = nexttile(t5); hold on; grid on;
        for ic = 1:nC
            e = pick(C, names{ic}, scales(is));
            if isempty(e) || isempty(e.traj) || all(isnan(e.traj.theta_hat(:)))
                continue;
            end
            est(ic) = true;
            plot(e.traj.t, vecnorm(e.traj.theta_hat,2,1), '-', ...
                 'Color', cols(ic,:), 'HandleVisibility', 'off');
            plot(e.traj.t, vecnorm(e.traj.theta_true,2,1), ':', ...
                 'Color', cols(ic,:), 'LineWidth', 1.0, 'HandleVisibility', 'off');
        end
        title(sprintf('model scale = %.2g', scales(is)));
    end
    % One legend row, colour for the controller and line style for estimate or
    % truth, as in fig 4: a row per (controller, quantity) pair stood two rows
    % tall over a column-width figure. LaTeX for the hat, so a controller
    % name's underscore is escaped (under 'tex' l1_con printed a subscript c).
    for ic = find(est)
        plot(ax5(1), NaN, NaN, '-', 'Color', cols(ic,:), 'LineWidth', 1.2, ...
             'DisplayName', strrep(names{ic}, '_', '\_'));
    end
    plot(ax5(1), NaN, NaN, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, ...
         'DisplayName', '$\|\hat{\theta}\|$');
    plot(ax5(1), NaN, NaN, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, ...
         'DisplayName', 'true $\|\theta\|$');
    ylabel(t5, '$\|\theta\|$', 'Interpreter', 'latex');
    xlabel(t5, 'Time (s)');
    % in a tile of its own: every footstrike's spike reaches the top of a panel
    lg = legend(ax5(1), 'Interpreter', 'latex', 'Orientation', 'horizontal');
    lg.Layout.Tile = 'north';
    share_time(ax5, [0 tmax]);
    figs(end+1) = f5;
end

% ------------------------------------------------------------------- saving
if ~isempty(savedir)
    if ~exist(savedir,'dir'), mkdir(savedir); end
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(savedir, sprintf('ch4_fig%d.png', k)), ...
                       'Resolution', 600);
    end
    fprintf(' ch4_plot_uncertainty: %d figures saved to %s\n', ...
            numel(figs), savedir);
end

end

% ---------------------------------------------------------------------------
function f = panel_grid(C, names, scales, cols, field, ylab, sc, ttl)
%PANEL_GRID  ny rows x nS columns of a per-output signal, cases across columns.
nS = numel(scales);
e0 = first_traj(C);
nr = size(e0.traj.(field), 1);

f = ch4_figure(ttl, 'page', 1.8 + 2.6*nr);
t = tiledlayout(f, nr, nS, 'TileSpacing', 'tight', 'Padding', 'compact');
ax = gobjects(nr, nS);
for is = 1:nS
    for ir = 1:nr
        ax(ir,is) = nexttile(t, (ir-1)*nS + is); hold on; grid on;
        for ic = 1:numel(names)
            e = pick(C, names{ic}, scales(is));
            if isempty(e) || isempty(e.traj), continue; end
            plot(e.traj.t, e.traj.(field)(ir,:) * sc, ...
                 'Color', cols(ic,:), 'DisplayName', names{ic});
        end
        if is == 1, ylabel(sprintf(ylab, ir)); end
        if ir == 1, title(sprintf('scale = %.2g', scales(is))); end
    end
end
lg = legend(ax(1,1), 'Interpreter', 'none', 'Orientation', 'horizontal');
lg.Layout.Tile = 'north';
% Share the y-axis ACROSS cases within each output row, so degradation with
% perturbation is visible rather than normalized away...
for ir = 1:nr, link_y(ax(ir,:)); end

% ...but then CLIP it. These are signed signals, so unlike the CLF panel they
% cannot go on a log axis, and a controller that has lost the robot swings
% +-300 deg while the ones that are working stay inside a few degrees. Autoscale
% shows the failure perfectly and every curve worth reading as a flat line on
% zero -- which is the wrong way round, since the failure is already legible
% from the table and the interesting question is what the survivors did.
%
% So the limits come from a high percentile of the row rather than its max, and
% the ylabel says so wherever a curve actually leaves the axes. Silent clipping
% would be worse than either extreme.
clipped = false;
for ir = 1:nr
    vals = [];
    for is = 1:nS
        for ic = 1:numel(names)
            e = pick(C, names{ic}, scales(is));
            if isempty(e) || isempty(e.traj), continue; end
            vals = [vals, e.traj.(field)(ir,:) * sc]; %#ok<AGROW>
        end
    end
    if isempty(vals), continue; end

    lim = prctile(abs(vals), 98) * 1.3;
    if ~isfinite(lim) || lim <= 0, continue; end

    if max(abs(vals)) > lim
        clipped = true;
        for is = 1:nS
            if isgraphics(ax(ir,is)), ylim(ax(ir,is), [-lim lim]); end
        end
        if isgraphics(ax(ir,1))
            ylabel(ax(ir,1), sprintf([ylab ' *'], ir));
        end
    end
end

% one footnote for the whole figure rather than per panel, under the shared
% time label, where no panel's decorations can cover it
note = {};
if clipped
    note = {['\fontsize{8}\color[rgb]{0.35,0.35,0.35}* axis clipped to the ' ...
             '98th percentile; a diverging run leaves the frame']};
end
xlabel(t, [{'Time (s)'}, note]);

% a column is one case, so its rows plot the same runs over the same time
for is = 1:nS, share_time(ax(:,is)); end
end

function e = pick(C, name, scale)
idx = find(strcmpi({C.name}, name) & [C.mass_scale] == scale, 1);
if isempty(idx), e = []; else, e = C(idx); end
end

function e = first_traj(C)
for k = 1:numel(C)
    if ~isempty(C(k).traj), e = C(k); return; end
end
error('ch4_plot_uncertainty:noTraj', 'No stored trajectories in the sweep.');
end

function link_y(ax)
ax = ax(isgraphics(ax));
if numel(ax) > 1, linkaxes(ax, 'y'); end
end

function share_time(ax, xl)
%SHARE_TIME  Stacked panels on one time axis, labelled on the bottom panel.
% The ticks are copied from the bottom panel rather than left to each axes: an
% axes whose tick labels are hidden picks its ticks differently, and its grid
% lines then disagreed with the panel below. Call it once the figure's labels
% and legends are in place, since the ticks MATLAB picks depend on the size
% the layout leaves the panel (read before drawnow, they came out as 0 and 5).
linkaxes(ax, 'x');
if nargin > 1, xlim(ax(1), xl); end
drawnow;
set(ax, 'XTick', ax(end).XTick);
set(ax(1:end-1), 'XTickLabel', {});
end

function decade_ticks(ax)
%DECADE_TICKS  Grid a log axis at every decade and label every other one.
% Left to itself, a panel this short keeps a single label, and a log axis with
% one label cannot be read.
e   = ceil(log10(ax.YLim(1))):floor(log10(ax.YLim(2)));
lab = compose('10^{%d}', e);
lab(mod(e, 2) ~= 0) = {''};
set(ax, 'YTick', 10.^e, 'YTickLabel', lab);
end

function s = roman(k)
r = {'I','II','III','IV','V','VI'};
if k <= numel(r), s = r{k}; else, s = num2str(k); end
end

function s = case_of_scale(s_mass)
%CASE_OF_SCALE  The chapter's case numeral for a mass scale.
%
% NUMBER THE PANEL BY ITS PERTURBATION, NOT BY ITS POSITION. Section 4.1.4
% numbers the scales I = 1, II = 1.5, III = 0.7, IV = 3, but the two presets
% sweep them in different orders -- 'robust' runs 1, 1.5, 0.7 and 'l1' runs
% 1, 0.7, 1.5 -- so numbering by position labelled the L1 panels II and III in
% the reverse of the robust table, and the same mass scale carried two
% different numerals across the report's figures.
switch round(100 * s_mass)
    case 100, s = 'I';
    case 150, s = 'II';
    case  70, s = 'III';
    case 300, s = 'IV';
    otherwise, s = sprintf('%.2g', s_mass);
end
end
