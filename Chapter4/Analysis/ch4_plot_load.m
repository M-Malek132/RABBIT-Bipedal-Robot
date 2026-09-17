function figs = ch4_plot_load(S, p, savedir, opts)
%CH4_PLOT_LOAD  The unknown-load figures (Fig. 4.11), from ch4_load_study.
%
%   figs = ch4_plot_load(S, p)
%   figs = ch4_plot_load(S, p, savedir)
%   figs = ch4_plot_load(S, p, savedir, opts)
%
%   fig 1  Fig. 4.11a, the random load. Left: one controller mid-step on each
%          of the first steps, with the mass it carries on that step drawn at
%          the hip. Top right: the load drawn for every step. Bottom right:
%          each controller's largest tracking error in every step -- a line
%          that ends early is a controller that fell.
%   fig 2  Fig. 4.11b, the norm of the joint torques under each fixed load, one
%          panel per controller, over the first steps.
%   fig 3  torso phase portrait of one controller under each fixed load, over
%          the nominal orbit. A uniform mass scale cannot move the orbit (see
%          ch4_plot_uncertainty); a load changes M non-uniformly and can, so
%          this is where the two perturbations part ways.
%
% Inputs
%   S       : struct array from ch4_load_study (needs .traj)
%   p       : parameter struct
%   savedir : optional directory; figures are written as ch4_fig*.png
%   opts    : struct with
%               .controller  whose snapshots and phase portrait (default 'l1_con')
%               .snap_steps  steps in the snapshot strip (default 5, as Fig. 4.11a)
%               .show_steps  steps in the torque panels (default 4)
%               .X_orbit     nominal-orbit states for fig 3 (optional)
%
% Output
%   figs : vector of figure handles
%
% See also CH4_LOAD_STUDY, CH4_DRAW_ROBOT, CH4_PLOT_UNCERTAINTY.

if nargin < 3, savedir = ''; end
if nargin < 4, opts = struct(); end
opts = fill_defaults(opts, struct('controller', 'l1_con', 'snap_steps', 5, ...
                                  'show_steps', 4, 'X_orbit', []));

if isempty(S) || all(cellfun(@isempty, {S.traj}))
    error('ch4_plot_load:noTraj', ...
          'The study carries no trajectories. Re-run with opts.store_traj = true.');
end

names   = unique({S.name}, 'stable');
nC      = numel(names);
cols    = lines(max(nC, 3));
is_rand = isnan([S.load_mass]);
figs    = gobjects(0);

% one colour per fixed load, light to heavy
load_cols = [0.30 0.60 0.90; 0.95 0.60 0.10; 0.75 0.15 0.15; 0.40 0.40 0.40];

%% --------------------------------------------------- fig 1: the random load
R = S(is_rand);
if ~isempty(R)
    % Snapshots beside the per-step panels, not above them: the robot is
    % taller than its stride, so a strip across the page drew five small
    % robots in a band of white.
    f1 = ch4_figure('ch4: random unknown load', 'page', 8.6);
    t1 = tiledlayout(f1, 2, 5, 'TileSpacing', 'tight', 'Padding', 'compact');

    ax1 = nexttile(t1, 1, [2 2]); hold(ax1, 'on'); axis(ax1, 'equal');
    e = pick(R, opts.controller);
    if ~isempty(e) && ~isempty(e.traj)
        tb = [0, cumsum(e.step_T)];
        n  = min(opts.snap_steps, numel(e.step_T));
        xs = [];
        for k = 1:n
            % mid-step: the legs pass each other, so neighbouring snapshots
            % overlap least
            [~, i] = min(abs(e.traj.t - 0.5 * (tb(k) + tb(k+1))));
            b  = ch4_draw_robot(ax1, e.traj.x(1:p.nq, i), false);
            mk = 6 + 1.6 * sqrt(e.loads(k));        % marker area grows with mass
            plot(ax1, b.hip(1), b.hip(2), 'o', 'MarkerSize', mk, ...
                 'MarkerFaceColor', [0.95 0.7 0.2], 'MarkerEdgeColor', [0.5 0.35 0]);
            text(ax1, b.torso_top(1), b.torso_top(2) + 0.12, ...
                 sprintf('%.0f kg', e.loads(k)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8);
            xs = [xs, b.stance_foot(1), b.swing_foot(1), b.hip(1)]; %#ok<AGROW>
        end
        xl = [min(xs) - 0.3, max(xs) + 0.3];
        plot(ax1, xl, [0 0], 'k-', 'LineWidth', 1.5);
        xlim(ax1, xl); ylim(ax1, [-0.05, 1.95]);
    end
    % the caption says what the mass is (never told to the controller, redrawn
    % every step); a title saying it was wider than the panel
    title(ax1, sprintf('%s at mid-step, steps 1-%d', ...
                       opts.controller, opts.snap_steps), 'Interpreter', 'none');
    xlabel(ax1, 'x (m)'); ylabel(ax1, 'z (m)');

    % the draws: identical for every controller (ch4_simulate seeds them), but
    % a run that fell stops drawing, so take the longest sequence
    [~, il] = max(arrayfun(@(r) numel(r.loads), R));
    loads = R(il).loads;
    ns    = numel(loads);

    ax2 = nexttile(t1, 3, [1 3]); hold(ax2, 'on'); grid(ax2, 'on');
    bar(ax2, 1:ns, loads, 0.6, 'FaceColor', [0.95 0.7 0.2], 'EdgeColor', 'none');
    xlim(ax2, [0.5, ns + 0.5]);
    % left alone, the axis kept two labels (0 and 50), too few to read a load by
    ytop = 10 * ceil(max(loads) / 10);
    ylim(ax2, [0, ytop]); yticks(ax2, 0:20:ytop);
    ylabel(ax2, 'load (kg)');
    title(ax2, 'mass carried on each step');

    ax3 = nexttile(t1, 8, [1 3]); hold(ax3, 'on'); grid(ax3, 'on');
    for ic = 1:nC
        e = pick(R, names{ic});
        if isempty(e) || isempty(e.traj), continue; end
        pk = per_step_max(e.traj.x, e.traj.eta_n, e.steps_completed);
        plot(ax3, 1:numel(pk), pk, '-o', 'Color', cols(ic,:), 'LineWidth', 1.1, ...
             'MarkerSize', 3, 'MarkerFaceColor', cols(ic,:), ...
             'DisplayName', sprintf('%s (%d steps)', names{ic}, e.steps_completed));
    end
    xlim(ax3, [0.5, ns + 0.5]);
    xlabel(ax3, 'step'); ylabel(ax3, 'max ||\eta|| in step');
    legend(ax3, 'Location', 'northoutside', 'NumColumns', 2, ...
           'Interpreter', 'none');
    figs(end+1) = f1;
end

%% ---------------------------------------- fig 2: torques under fixed loads
Fx = S(~is_rand);
if ~isempty(Fx)
    loads = unique([Fx.load_mass], 'stable');

    f2 = ch4_figure('ch4: norm of torques under a fixed unknown load', 'page', 6.2);
    % 'compact', not 'tight': with the inner scales unlabelled, 'tight' butted
    % the panels together
    t2 = tiledlayout(f2, 1, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = gobjects(1, nC);
    vals = [];
    for ic = 1:nC
        ax(ic) = nexttile(t2); hold(ax(ic), 'on'); grid(ax(ic), 'on');
        for il = 1:numel(loads)
            e = pick(Fx, names{ic}, loads(il));
            if isempty(e) || isempty(e.traj), continue; end
            tb   = [0, cumsum(e.step_T)];
            tend = tb(min(opts.show_steps, numel(e.step_T)) + 1);
            in   = e.traj.t <= tend + 1e-9;
            un   = vecnorm(e.traj.u(:, in), 2, 1);
            vals = [vals, un]; %#ok<AGROW>
            plot(ax(ic), e.traj.t(in), un, 'Color', ...
                 load_cols(1 + mod(il - 1, size(load_cols, 1)), :), ...
                 'DisplayName', sprintf('%g kg', loads(il)));
        end
        title(ax(ic), names{ic}, 'Interpreter', 'none');
    end
    xlabel(t2, 'Time (s)'); ylabel(t2, 'norm of torques (Nm)');
    % inside the last panel's top right, clear of an adaptive law's torques; in
    % a tile of its own it took a centimetre from the panels
    legend(ax(end), 'Location', 'northeast');
    % Shared axes, so the panels compare -- clipped at a high percentile, as in
    % ch4_plot_uncertainty, so one controller's spike does not flatten the rest.
    ax = ax(isgraphics(ax));
    if numel(ax) > 1, linkaxes(ax, 'xy'); end
    if ~isempty(vals)
        lim = prctile(vals, 99) * 1.2;
        if isfinite(lim) && lim > 0
            ylim(ax(1), [0, lim]);
            if max(vals) > lim
                xlabel(t2, {'Time (s)', ['\fontsize{8}\color[rgb]{0.35,0.35,0.35}' ...
                       'axis clipped at 1.2 x the 99th percentile; spikes leave the frame']});
            end
        end
    end
    % The scale is labelled on the left panel only, and its ticks are copied to
    % the others: an axes whose tick labels are hidden picks its ticks
    % differently, and its grid lines then disagreed with the labelled panel.
    % drawnow first, because the ticks MATLAB picks depend on the panel's final
    % size in the layout.
    drawnow;
    set(ax, 'YTick', ax(1).YTick);
    set(ax(2:end), 'YTickLabel', {});
    figs(end+1) = f2;

    %% ------------------------------------ fig 3: does a load move the orbit?
    % Two panels on shared axes: the unconstrained baseline and the adaptive law.
    % What separates them at the same load is tracking error; what both share is
    % the load's effect on the robot. Not the boxed baseline: it falls within a
    % few steps under every load (see ch4_load_study), so its panel would show
    % falls rather than orbits, and its fall would set the axes for both.
    panels = unique({'clfqp', opts.controller}, 'stable');

    f3 = ch4_figure('ch4: torso phase portrait under a fixed load', 'page', 8.0);
    t3 = tiledlayout(f3, 1, numel(panels), 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = gobjects(1, numel(panels));
    for ip = 1:numel(panels)
        ax(ip) = nexttile(t3); hold(ax(ip), 'on'); grid(ax(ip), 'on');
        if ~isempty(opts.X_orbit)
            plot(ax(ip), opts.X_orbit(3, :) * 180/pi, opts.X_orbit(3 + p.nq, :) * 180/pi, ...
                 'k-', 'LineWidth', 1.4, 'DisplayName', 'nominal orbit, no load');
        end
        for il = 1:numel(loads)
            e = pick(Fx, panels{ip}, loads(il));
            if isempty(e) || isempty(e.traj), continue; end
            plot(ax(ip), e.traj.x(3, :) * 180/pi, e.traj.x(3 + p.nq, :) * 180/pi, '-', ...
                 'Color', load_cols(1 + mod(il - 1, size(load_cols, 1)), :), ...
                 'LineWidth', 0.6, 'DisplayName', sprintf('%g kg', loads(il)));
        end
        xlabel(ax(ip), 'q_{torso} (deg)');
        if ip == 1, ylabel(ax(ip), 'dq_{torso} (deg/s)'); end
        title(ax(ip), panels{ip}, 'Interpreter', 'none');
    end
    if numel(ax) > 1, linkaxes(ax, 'xy'); end
    % What the comparison means -- a load changes M non-uniformly, so unlike a
    % mass scale it can move the orbit, and the panels differ by tracking
    % error -- is the report caption's to say; as a subtitle it was a second
    % line of small print.
    title(t3, 'Torso phase portrait under a fixed unknown load', ...
          'FontSize', 9.9, 'FontWeight', 'bold');
    % in the baseline panel's empty top, as in fig 2
    legend(ax(1), 'Interpreter', 'none', 'NumColumns', 2, 'Location', 'north');
    drawnow;                                  % ticks as in fig 2
    set(ax, 'YTick', ax(1).YTick);
    set(ax(2:end), 'YTickLabel', {});
    figs(end+1) = f3;
end

%% ------------------------------------------------------------------- saving
if ~isempty(savedir)
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(savedir, sprintf('ch4_fig%d.png', k)), ...
                       'Resolution', 600);
    end
    fprintf(' ch4_plot_load: %d figures saved to %s\n', numel(figs), savedir);
end

end

% ---------------------------------------------------------------------------
function e = pick(S, name, load_mass)
sel = strcmpi({S.name}, name);
if nargin > 2, sel = sel & [S.load_mass] == load_mass; end
idx = find(sel, 1);
if isempty(idx), e = []; else, e = S(idx); end
end

function pk = per_step_max(X, v, n_steps)
%PER_STEP_MAX  Largest v within each completed step.
%
% SPLIT ON THE RELABEL, NOT ON TIME. Every footstrike puts two samples at one
% instant, the pre-impact state closing a step and the post-impact state
% opening the next, and decimation may keep either or both. A time window
% therefore hands the post-impact sample -- usually the step's largest error,
% since the impact is what throws the outputs off -- to the step that ended,
% and neighbouring steps read the same maximum. At impact the legs swap roles,
% so the stance-hip angle jumps by the gap between the legs; that jump is
% where a step begins, whichever samples survived.
starts = [1, find(abs(diff(X(4, :))) > 0.2) + 1, size(X, 2) + 1];
pk = nan(1, n_steps);
if numel(starts) - 1 ~= n_steps
    warning('ch4_plot_load:steps', ...
            'Found %d relabels for %d steps; per-step maxima omitted.', ...
            numel(starts) - 2, n_steps);
    return;
end
for k = 1:n_steps
    pk(k) = max(v(starts(k):starts(k+1) - 1));
end
end

function s = fill_defaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i})
        s.(f{i}) = d.(f{i});
    end
end
end
