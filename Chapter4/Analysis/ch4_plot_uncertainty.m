function figs = ch4_plot_uncertainty(C, p, savedir)
%CH4_PLOT_UNCERTAINTY  The Chapter-4 figures, from a comparison sweep.
%
%   figs = ch4_plot_uncertainty(C, p)
%   figs = ch4_plot_uncertainty(C, p, savedir)
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
%   p       : parameter struct, for labels
%   savedir : optional directory; figures are written as ch4_fig*.png
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

cols = lines(max(nC,3));
figs = gobjects(0);

% ---------------------------------------------------------------- fig 1: CLF
f1 = figure('Name','ch4: CLF under model perturbation', ...
            'Position',[80 80 620 240*nS]);
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
    ax(is) = subplot(nS,1,is); hold on; grid on;
    for ic = 1:nC
        e = pick(C, names{ic}, scales(is));
        if isempty(e) || isempty(e.traj), continue; end
        Vv = e.traj.V;
        if all(isnan(Vv)), Vv = e.traj.eta_n.^2; end
        plot(e.traj.t, max(Vv, Vfloor), 'Color', cols(ic,:), ...
             'LineWidth', 1.3, 'DisplayName', names{ic});
    end
    set(gca, 'YScale', 'log');
    ylim([Vfloor, Vmax*2]);
    ylabel('V_\epsilon');
    title(sprintf('Case %s: model scale = %.2g', roman(is), scales(is)));
    if is == 1, legend('Location','southeast','Interpreter','none'); end
    if is == nS, xlabel('Time (s)'); end
end
figs(end+1) = f1;

% ------------------------------------------------- fig 2 / 3: outputs, torques
f2 = panel_grid(C, names, scales, cols, 'y',  'y_%d (deg)', 180/pi, ...
                'ch4: tracking errors');
f3 = panel_grid(C, names, scales, cols, 'u',  'u_%d (Nm)',  1, ...
                'ch4: joint torques');
figs(end+1) = f2;
figs(end+1) = f3;

% ------------------------------------------------------ fig 4: phase portrait
f4 = figure('Name','ch4: torso phase portrait','Position',[120 120 640 420]);
hold on; grid on;
mk = {'-','--',':','-.'};
for ic = 1:nC
    for is = 1:nS
        e = pick(C, names{ic}, scales(is));
        if isempty(e) || isempty(e.traj), continue; end
        qt  = e.traj.x(3,  :) * 180/pi;
        dqt = e.traj.x(3+p.nq, :) * 180/pi;
        plot(qt, dqt, mk{min(is,numel(mk))}, 'Color', cols(ic,:), ...
             'LineWidth', 1.1, ...
             'DisplayName', sprintf('%s, scale %.2g', names{ic}, scales(is)));
    end
end
xlabel('q_{torso} (deg)'); ylabel('dq_{torso} (deg/s)');
title('Torso phase portrait: uncertainty moves the periodic orbit');
legend('Location','best','Interpreter','none');
figs(end+1) = f4;

% ------------------------------------------------------- fig 5: L1 estimator
has_l1 = any(cellfun(@(n) any(strcmpi(n,{'l1','l1_con'})), names));
if has_l1
    f5 = figure('Name','ch4: L1 estimator','Position',[160 160 620 240*nS]);
    for is = 1:nS
        subplot(nS,1,is); hold on; grid on;
        for ic = 1:nC
            e = pick(C, names{ic}, scales(is));
            if isempty(e) || isempty(e.traj) || all(isnan(e.traj.theta_hat(:)))
                continue;
            end
            plot(e.traj.t, vecnorm(e.traj.theta_hat,2,1), '-', ...
                 'Color', cols(ic,:), 'LineWidth', 1.3, ...
                 'DisplayName', sprintf('%s: ||\\theta hat||', names{ic}));
            plot(e.traj.t, vecnorm(e.traj.theta_true,2,1), ':', ...
                 'Color', cols(ic,:), 'LineWidth', 1.1, ...
                 'DisplayName', sprintf('%s: ||\\theta true||', names{ic}));
        end
        ylabel('||\theta||');
        title(sprintf('model scale = %.2g', scales(is)));
        if is == 1, legend('Location','northeast','Interpreter','tex'); end
        if is == nS, xlabel('Time (s)'); end
    end
    figs(end+1) = f5;
end

% ------------------------------------------------------------------- saving
if ~isempty(savedir)
    if ~exist(savedir,'dir'), mkdir(savedir); end
    for k = 1:numel(figs)
        saveas(figs(k), fullfile(savedir, sprintf('ch4_fig%d.png', k)));
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

f = figure('Name', ttl, 'Position', [100 100 320*nS 170*nr]);
ax = gobjects(nr, nS);
for is = 1:nS
    for ir = 1:nr
        ax(ir,is) = subplot(nr, nS, (ir-1)*nS + is); hold on; grid on;
        for ic = 1:numel(names)
            e = pick(C, names{ic}, scales(is));
            if isempty(e) || isempty(e.traj), continue; end
            plot(e.traj.t, e.traj.(field)(ir,:) * sc, ...
                 'Color', cols(ic,:), 'LineWidth', 1.1, ...
                 'DisplayName', names{ic});
        end
        if is == 1, ylabel(sprintf(ylab, ir)); end
        if ir == 1
            title(sprintf('scale = %.2g', scales(is)));
            if is == nS
                legend('Location','best','Interpreter','none','FontSize',7);
            end
        end
        if ir == nr, xlabel('Time (s)'); end
    end
end
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
        for is = 1:nS
            if isgraphics(ax(ir,is)), ylim(ax(ir,is), [-lim lim]); end
        end
        if isgraphics(ax(ir,1))
            ylabel(ax(ir,1), sprintf([ylab ' *'], ir));
        end
    end
end

% one footnote for the whole figure rather than per panel
annotation(f, 'textbox', [0.005 0.005 0.5 0.03], 'String', ...
    '* axis clipped to the 98th percentile; a diverging run leaves the frame', ...
    'EdgeColor', 'none', 'FontSize', 7, 'Color', [0.35 0.35 0.35]);
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

function s = roman(k)
r = {'I','II','III','IV','V','VI'};
if k <= numel(r), s = r{k}; else, s = num2str(k); end
end
