function figs = ch4_plot_load(S, p, savedir, opts)
%CH4_PLOT_LOAD  The unknown-load figures (Fig. 4.11), from ch4_load_study.
%
%   figs = ch4_plot_load(S, p)
%   figs = ch4_plot_load(S, p, savedir)
%   figs = ch4_plot_load(S, p, savedir, opts)
%
%   fig 1  Fig. 4.11a, the random load. Top: one controller mid-step on each of
%          the first steps, with the mass it carries on that step drawn at the
%          hip. Middle: the load drawn for every step. Bottom: each controller's
%          largest tracking error in every step -- a line that ends early is a
%          controller that fell.
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
    f1 = figure('Name', 'ch4: random unknown load', 'Position', [60 60 1000 860]);

    ax1 = subplot(3, 1, 1); hold(ax1, 'on'); axis(ax1, 'equal');
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
                 'HorizontalAlignment', 'center', 'FontSize', 9);
            xs = [xs, b.stance_foot(1), b.swing_foot(1), b.hip(1)]; %#ok<AGROW>
        end
        xl = [min(xs) - 0.3, max(xs) + 0.3];
        plot(ax1, xl, [0 0], 'k-', 'LineWidth', 1.5);
        xlim(ax1, xl); ylim(ax1, [-0.05, 1.95]);
    end
    title(ax1, sprintf(['%s carrying a mass it is never told about, redrawn ' ...
                        'every step: steps 1-%d, mid-step'], ...
                       opts.controller, opts.snap_steps), 'Interpreter', 'none');
    xlabel(ax1, 'x [m]'); ylabel(ax1, 'z [m]');

    % the draws: identical for every controller (ch4_simulate seeds them), but
    % a run that fell stops drawing, so take the longest sequence
    [~, il] = max(arrayfun(@(r) numel(r.loads), R));
    loads = R(il).loads;
    ns    = numel(loads);

    ax2 = subplot(3, 1, 2); hold(ax2, 'on'); grid(ax2, 'on');
    bar(ax2, 1:ns, loads, 0.6, 'FaceColor', [0.95 0.7 0.2], 'EdgeColor', 'none');
    xlim(ax2, [0.5, ns + 0.5]);
    ylabel(ax2, 'load [kg]');
    title(ax2, 'mass carried on each step');

    ax3 = subplot(3, 1, 3); hold(ax3, 'on'); grid(ax3, 'on');
    for ic = 1:nC
        e = pick(R, names{ic});
        if isempty(e) || isempty(e.traj), continue; end
        pk = per_step_max(e.traj.x, e.traj.eta_n, e.steps_completed);
        plot(ax3, 1:numel(pk), pk, '-o', 'Color', cols(ic,:), 'LineWidth', 1.3, ...
             'MarkerSize', 4, 'MarkerFaceColor', cols(ic,:), ...
             'DisplayName', sprintf('%s (%d steps)', names{ic}, e.steps_completed));
    end
    xlim(ax3, [0.5, ns + 0.5]);
    xlabel(ax3, 'step'); ylabel(ax3, 'max ||\eta|| in step');
    legend(ax3, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
           'Interpreter', 'none');
    figs(end+1) = f1;
end

%% ---------------------------------------- fig 2: torques under fixed loads
Fx = S(~is_rand);
if ~isempty(Fx)
    loads = unique([Fx.load_mass], 'stable');

    f2 = figure('Name', 'ch4: norm of torques under a fixed unknown load', ...
                'Position', [80 80 400*nC 380]);
    ax = gobjects(1, nC);
    vals = [];
    for ic = 1:nC
        ax(ic) = subplot(1, nC, ic); hold(ax(ic), 'on'); grid(ax(ic), 'on');
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
                 'LineWidth', 1.2, 'DisplayName', sprintf('%g kg', loads(il)));
        end
        title(ax(ic), names{ic}, 'Interpreter', 'none');
        xlabel(ax(ic), 'Time (s)');
        if ic == 1, ylabel(ax(ic), 'norm of torques (Nm)'); end
        if ic == nC, legend(ax(ic), 'Location', 'best'); end
    end
    % Shared axis, so the panels compare -- clipped at a high percentile, as in
    % ch4_plot_uncertainty, so one controller's spike does not flatten the rest.
    ax = ax(isgraphics(ax));
    if numel(ax) > 1, linkaxes(ax, 'y'); end
    if ~isempty(vals)
        lim = prctile(vals, 99) * 1.2;
        if isfinite(lim) && lim > 0
            ylim(ax(1), [0, lim]);
            if max(vals) > lim
                annotation(f2, 'textbox', [0.005 0.005 0.6 0.04], 'String', ...
                    'axis clipped at 1.2 x the 99th percentile; spikes leave the frame', ...
                    'EdgeColor', 'none', 'FontSize', 7, 'Color', [0.35 0.35 0.35]);
            end
        end
    end
    figs(end+1) = f2;

    %% ------------------------------------ fig 3: does a load move the orbit?
    % Two panels on shared axes: the unconstrained baseline and the adaptive law.
    % What separates them at the same load is tracking error; what both share is
    % the load's effect on the robot. Not the boxed baseline: it falls within a
    % few steps under every load (see ch4_load_study), so its panel would show
    % falls rather than orbits, and its fall would set the axes for both.
    panels = unique({'clfqp', opts.controller}, 'stable');

    f3 = figure('Name', 'ch4: torso phase portrait under a fixed load', ...
                'Position', [120 120 560*numel(panels) 460]);
    ax = gobjects(1, numel(panels));
    for ip = 1:numel(panels)
        ax(ip) = subplot(1, numel(panels), ip); hold(ax(ip), 'on'); grid(ax(ip), 'on');
        if ~isempty(opts.X_orbit)
            plot(ax(ip), opts.X_orbit(3, :) * 180/pi, opts.X_orbit(3 + p.nq, :) * 180/pi, ...
                 'k-', 'LineWidth', 2, 'DisplayName', 'nominal orbit, no load');
        end
        for il = 1:numel(loads)
            e = pick(Fx, panels{ip}, loads(il));
            if isempty(e) || isempty(e.traj), continue; end
            plot(ax(ip), e.traj.x(3, :) * 180/pi, e.traj.x(3 + p.nq, :) * 180/pi, '-', ...
                 'Color', load_cols(1 + mod(il - 1, size(load_cols, 1)), :), ...
                 'LineWidth', 1, 'DisplayName', sprintf('%g kg', loads(il)));
        end
        xlabel(ax(ip), 'q_{torso} (deg)');
        if ip == 1, ylabel(ax(ip), 'dq_{torso} (deg/s)'); end
        title(ax(ip), panels{ip}, 'Interpreter', 'none');
        legend(ax(ip), 'Location', 'best', 'Interpreter', 'none');
    end
    if numel(ax) > 1, linkaxes(ax, 'xy'); end
    sgtitle(f3, {'Torso phase portrait under a fixed unknown load', ...
                 ['a load changes M non-uniformly, so unlike a mass scale it can ' ...
                  'move the orbit; the two panels differ by tracking error']});
    figs(end+1) = f3;
end

%% ------------------------------------------------------------------- saving
if ~isempty(savedir)
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    for k = 1:numel(figs)
        saveas(figs(k), fullfile(savedir, sprintf('ch4_fig%d.png', k)));
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
