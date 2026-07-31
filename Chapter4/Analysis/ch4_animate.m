function out = ch4_animate(x0, alpha, p, controllers, n_steps, gifpath)
%CH4_ANIMATE  Animate the SAME perturbed robot under different controllers.
%
%   ch4_animate(x0, alpha, p)
%   ch4_animate(x0, alpha, p, {'clfqp','rclfqp_con'}, 4)
%   out = ch4_animate(x0, alpha, p, ctrls, 4, 'Results/ch4_walk.gif')
%
% ch3_animate shows one gait running. The Chapter-4 question is different and
% needs a different picture: given ONE robot that is not the model, which
% controllers keep it walking? So this draws each controller in its own panel,
% on the same perturbed plant, SYNCHRONIZED IN TIME.
%
% WHY SYNCHRONIZED, AND WHY THAT TAKES WORK.  The runs do not share a time grid
% -- each has its own adaptive solver output -- and, more importantly, they do
% not share a DURATION, because a controller that loses the robot stops early.
% Playing each at its own pace would put the panels out of step and make the
% comparison meaningless just when it gets interesting. So every run is
% interpolated onto one common grid, and a run that has already ended is frozen
% at its last pose and labelled. The frozen panel IS the result: that controller
% is on the floor while the others are still walking.
%
% Interpolation is on q only, and only for drawing. Nothing here feeds back into
% a simulation, so it cannot launder a bad trajectory into a good-looking one --
% and a run that ended early is marked rather than extrapolated.
%
% Each panel reports its live ||eta||, because two stick figures can look
% similar for a while before one of them falls, and the tracking error is what
% distinguishes them well before that.
%
% Inputs
%   x0, alpha   : the gait (ch4_load_gait)
%   p           : parameter struct; p.uncertainty is what all runs share
%   controllers : cell of controller names. Default {'clfqp','rclfqp_con','l1'}
%   n_steps     : steps to attempt (default 4)
%   gifpath     : optional .gif to write
%
% Output
%   out : struct array, one per controller, with .name .sim .fell .T
%
% Requires a normal (JVM-enabled) MATLAB session for the figure.
%
% See also CH3_ANIMATE, CH4_SIMULATE, CH4_COMPARE_CONTROLLERS, CH3_BODY_POINTS.

if nargin < 4 || isempty(controllers)
    controllers = {'clfqp', 'rclfqp_con', 'l1'};
end
if nargin < 5 || isempty(n_steps), n_steps = 4; end
if nargin < 6, gifpath = ''; end

nC = numel(controllers);

%% --- run them all on the same perturbed plant ---------------------------
% Field list must match the assignment below EXACTLY -- MATLAB rejects
% runs(ic) = struct(...) with "dissimilar structures" if it does not, and the
% error names neither the missing field nor the line that diverged.
runs = struct('name', {}, 'sim', {}, 'fell', {}, 'T', {}, 'eta', {}, ...
              't_u', {}, 'i_u', {});

fprintf('ch4_animate: mass scale %.2f, load %.1f kg\n', ...
        p.uncertainty.mass_scale, p.uncertainty.load_mass);

for ic = 1:nC
    pc = p;
    pc.controller = controllers{ic};
    % Only the laws whose formulation carries a box are told about one, exactly
    % as in ch4_compare_controllers -- otherwise the picture would not match
    % the table.
    pc.limits.enable.torque = any(strcmpi(controllers{ic}, ...
                                          {'clfqp_con', 'rclfqp_con'}));

    sim = ch4_simulate(x0, alpha, pc, n_steps);

    % tracking error along the run, for the panel readout
    eta_n = nan(1, size(sim.x, 2));
    for k = 1:size(sim.x, 2)
        [yk, ydk] = ch3_outputs(sim.x(:,k), alpha, pc);
        eta_n(k)  = norm([yk; ydk]);
    end

    % NOT a ternary helper. MATLAB evaluates every argument BEFORE entering the
    % function, so ternary(isempty(t), 0, t(end)) indexes an empty array and
    % errors on exactly the case it was written to handle -- a controller that
    % completed no steps, which at mass scale 0.7 is the baseline and is the
    % whole point of the picture.
    if isempty(sim.t)
        T_run = 0;
    else
        T_run = sim.t(end);
    end

    % One sorted-unique time base per run, built once. sim.t repeats values at
    % every step boundary (each step starts at the previous step's end time),
    % and interp1 rejects non-unique sample points.
    [t_u, i_u] = unique(sim.t(:).');

    runs(ic) = struct('name', controllers{ic}, 'sim', sim, ...
                      'fell', sim.n_ok < n_steps, 'T', T_run, ...
                      'eta', eta_n, 't_u', t_u, 'i_u', i_u); %#ok<AGROW>

    if runs(ic).fell, note = '   <- stopped early'; else, note = ''; end
    fprintf('  %-11s %d/%d steps, %.3f s%s\n', controllers{ic}, sim.n_ok, ...
            n_steps, T_run, note);
end

if ~any([runs.T] > 0)
    error('ch4_animate:noSteps', 'No controller completed any step.');
end

%% --- one common clock ---------------------------------------------------
T_end  = max([runs.T]);
% 150 frames at the 0.04 s delay written below is a ~6 s loop -- smooth enough
% to read the gait, and about 40% smaller on disk than 240. These GIFs get
% committed, so the frame count is a repo-size decision as much as a visual one.
nframe = 150;
tg     = linspace(0, T_end, nframe);

%% --- figure -------------------------------------------------------------
allx = [];
for ic = 1:nC
    if ~isempty(runs(ic).sim.x), allx = [allx, runs(ic).sim.x(1,:)]; end %#ok<AGROW>
end

fig = figure('Color','w','Position',[80 80 480*nC 460]);
ax  = gobjects(1, nC);
for ic = 1:nC
    ax(ic) = subplot(1, nC, ic);
    hold(ax(ic),'on'); grid(ax(ic),'on'); axis(ax(ic),'equal');
    xlim(ax(ic), [min(allx)-0.8, max(allx)+0.8]);
    ylim(ax(ic), [-0.15, 1.6]);
    xlabel(ax(ic),'x [m]');
    if ic == 1, ylabel(ax(ic),'z [m]'); end
end

sg = sgtitle(fig, '');

first = true;
for f = 1:nframe
    for ic = 1:nC
        r = runs(ic);
        cla(ax(ic));
        plot(ax(ic), xlim(ax(ic)), [0 0], 'k-', 'LineWidth', 2);

        if isempty(r.sim.x)
            title(ax(ic), sprintf('%s  --  NEVER COMPLETED A STEP', r.name), ...
                  'Interpreter','none', 'Color', [0.75 0.1 0.1]);
            continue;
        end

        % freeze at the last pose once this run has ended
        done = tg(f) >= r.T;
        tq   = min(tg(f), r.T);
        q    = sample_at(r.t_u, r.sim.x(1:p.nq, r.i_u), tq);
        e    = sample_at(r.t_u, r.eta(r.i_u),           tq);

        draw_robot(ax(ic), q, done);

        if done && r.fell
            lbl = sprintf('%s  --  STOPPED at %.2f s', r.name, r.T);
            col = [0.75 0.1 0.1];
        else
            lbl = sprintf('%s  |  ||eta|| = %.3f', r.name, e);
            col = [0 0 0];
        end
        title(ax(ic), lbl, 'Interpreter','none', 'Color', col);
    end

    sg.String = sprintf(['mass scale %.2f, load %.1f kg   |   t = %.2f s   ' ...
                         '|   controller sees the NOMINAL model'], ...
                        p.uncertainty.mass_scale, p.uncertainty.load_mass, tg(f));
    drawnow limitrate;

    if ~isempty(gifpath)
        fr = getframe(fig);
        [ind, cm] = rgb2ind(frame2im(fr), 256);
        if first
            imwrite(ind, cm, gifpath, 'gif', 'LoopCount', Inf, 'DelayTime', 0.04);
            first = false;
        else
            imwrite(ind, cm, gifpath, 'gif', 'WriteMode','append', 'DelayTime', 0.04);
        end
    end
end

if ~isempty(gifpath)
    fprintf('ch4_animate: wrote %s\n', gifpath);
end

if nargout > 0
    out = rmfield(runs, {'eta', 't_u', 'i_u'});   % drawing scratch, not results
end

end

% ---------------------------------------------------------------------------
function draw_robot(ax, q, dimmed)
%DRAW_ROBOT  The stick figure, in ch3_animate's colours.
%
% A finished run is drawn washed out rather than removed: the pose still carries
% information (where it ended up, and in what attitude), but it must not read as
% a robot that is still walking.
b = ch3_body_points(q);

if dimmed
    c_st = [0.85 0.6 0.6]; c_sw = [0.6 0.68 0.85]; c_k = [0.55 0.55 0.55];
    lw = 1.5;
else
    c_st = [0.85 0.2 0.2]; c_sw = [0.2 0.4 0.85]; c_k = [0 0 0];
    lw = 2.5;
end

plot(ax, [b.hip(1) b.torso_top(1)], [b.hip(2) b.torso_top(2)], '-', ...
     'Color', c_k, 'LineWidth', lw+0.5);
plot(ax, [b.hip(1) b.stance_knee(1) b.stance_foot(1)], ...
         [b.hip(2) b.stance_knee(2) b.stance_foot(2)], '-o', ...
     'Color', c_st, 'LineWidth', lw, 'MarkerSize', 5);
plot(ax, [b.hip(1) b.swing_knee(1) b.swing_foot(1)], ...
         [b.hip(2) b.swing_knee(2) b.swing_foot(2)], '--s', ...
     'Color', c_sw, 'LineWidth', lw-0.5, 'MarkerSize', 5);
plot(ax, b.hip(1), b.hip(2), 'o', 'Color', c_k, ...
     'MarkerFaceColor', c_k, 'MarkerSize', 6);
end

function v = sample_at(t_u, Y, tq)
%SAMPLE_AT  Rows of Y at time tq, on an already sorted-unique time base.
%
% Y is either nq x nt (a pose) or 1 x nt (a scalar signal); both come back as a
% column. The dedup happened once per run in the caller rather than per frame --
% 240 frames x 3 controllers x 2 scales is 1440 calls, and `unique` on a 4000
% point vector is not free.
%
% Clamped, never extrapolated. tq is already limited to the run's own duration,
% so 'extrap' would only ever fire on floating-point overshoot at the endpoint,
% and silently inventing a pose past the end of a run is exactly the thing this
% animation must not do.
if isvector(Y), Y = Y(:).'; end

if numel(t_u) < 2
    v = Y(:,1);
    return;
end

tq = min(max(tq, t_u(1)), t_u(end));
v  = interp1(t_u, Y.', tq, 'linear').';
end
