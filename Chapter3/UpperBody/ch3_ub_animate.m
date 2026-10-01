function ch3_ub_animate(gait_file, n_steps, gifpath, speed)
%CH3_UB_ANIMATE  Animate a gait solved with the upper-body (crutch) wrench.
%a
%   ch3_ub_animate('Results/ch3_ub_v100_w500.mat')
%   ch3_ub_animate(file, n_steps)
%   ch3_ub_animate(file, n_steps, 'Results/crutch_w500.gif')
%   ch3_ub_animate(file, n_steps, gifpath, speed)     speed 1 = real time
%
% Simulates the saved gait WITH its wrench (p.ub) under the saved controller
% and draws, per frame:
%
%   * the stick figure (stance leg solid red, swing leg dashed blue);
%   * the CRUTCH: a grey shaft from the shoulder point to the ground, laid
%     along the force direction -- a crutch can only push along its shaft,
%     so the shaft is where the force says it is. Hidden when |F| < 5 N;
%   * the force itself as a green arrow at the shoulder (100 N = 0.25 m);
%   * a lower panel with Fx and Fz over time and a cursor.
%
% This is SIMULATION, not the collocation nodes, so a gait that only exists
% on the mesh falls over here.
%
% Accepts files saved as {z, o} (sweeps) or {zC3, oC3} / {zC2, oC2}.
%
% See also CH3_ANIMATE, CH3_UB_WRENCH, CH3_SIMULATE.

if nargin < 2 || isempty(n_steps), n_steps = 6;  end
if nargin < 3, gifpath = ''; end
if nargin < 4 || isempty(speed),   speed = 1;    end

[z, p] = load_gait(gait_file);
[X, ~, alpha] = ch3_col_unpack(z, p);

sim = ch3_simulate(X(:, 1), alpha, p, n_steps);
if sim.n_ok == 0
    error('ch3_ub_animate:noSteps', 'No steps completed: %s', sim.reason);
end
fprintf('ch3_ub_animate: %d/%d steps (%s)\n', sim.n_ok, n_steps, sim.reason);

% ---- frames on a uniform clock, never interpolated across an impact -------
DT    = 0.02;                                   % s of simulated time per frame
t_off = cumsum([0, arrayfun(@(s) s.T, sim.steps)]);
t_fr  = 0:DT:t_off(end);
nf    = numel(t_fr);
Xf    = zeros(p.nx, nf);
for i = 1:nf
    k = find(t_fr(i) >= t_off(1:end-1), 1, 'last');
    s = sim.steps(k);
    Xf(:, i) = interp1(s.t(:), s.x.', min(t_fr(i) - t_off(k), s.T)).';
end

% ---- wrench and shoulder point at every frame ------------------------------
W  = zeros(3, nf);  Rsh = zeros(2, nf);
for i = 1:nf
    s_ph = ch3_phase(Xf(:, i), p);
    [~, W(:, i), Rsh(:, i)] = ch3_ub_wrench(Xf(1:p.nq, i), s_ph, p);
end

% ---- figure ----------------------------------------------------------------
top = 0;
for i = 1:nf
    b = ch3_body_points(Xf(1:p.nq, i));  top = max(top, b.torso_top(2));
end
xs  = Xf(1, :);
fig = figure('Color', 'w', 'Position', [100 100 950 620]);
ax  = subplot(3, 1, [1 2]);  hold(ax, 'on');  axis(ax, 'equal');  grid(ax, 'on');
xlim(ax, [min(xs) - 0.8, max(xs) + 1.0]);  ylim(ax, [-0.1, top + 0.1]);
xlabel(ax, 'x [m]');  ylabel(ax, 'z [m]');

ax2 = subplot(3, 1, 3);  hold(ax2, 'on');  grid(ax2, 'on');
plot(ax2, t_fr, W(1, :), 'LineWidth', 1.5);
plot(ax2, t_fr, W(2, :), 'LineWidth', 1.5);
for k = 2:numel(t_off) - 1, xline(ax2, t_off(k), ':', 'Color', [0.6 0.6 0.6]); end
cur = xline(ax2, 0, 'k-', 'LineWidth', 1.2);
xlim(ax2, [0 t_off(end)]);  xlabel(ax2, 't [s]');  ylabel(ax2, 'crutch force [N]');
legend(ax2, {'F_x', 'F_z'}, 'Location', 'eastoutside');

E  = ch3_col_eval(z, p);  [~, T] = ch3_col_unpack(z, p);
ttl = sprintf('%s   |   v = %.2f m/s,  step %.2f m   |   w = %g', ...
              strrep(gait_file, '_', '\_'), E.L_step / T, E.L_step, p.ub.w);

A_SCALE = 0.25 / 100;                          % m per N for the force arrow
first = true;
for k = 1:nf
    cla(ax);
    plot(ax, xlim(ax), [0 0], 'k-', 'LineWidth', 2);
    b  = ch3_body_points(Xf(1:p.nq, k));
    F  = W(1:2, k);  sh = Rsh(:, k);

    % crutch shaft: shoulder -> ground along the force direction
    if norm(F) > 5 && F(2) > 1e-6
        d   = -F / norm(F);                    % from shoulder toward the tip
        tip = sh + d * (sh(2) / -d(2));        % where the shaft meets z = 0
        plot(ax, [sh(1) tip(1)], [sh(2) tip(2)], '-', 'Color', [0.55 0.55 0.55], ...
             'LineWidth', 4);
        plot(ax, tip(1), tip(2), 'o', 'Color', [0.4 0.4 0.4], ...
             'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6);
    end

    plot(ax, [b.hip(1) b.torso_top(1)], [b.hip(2) b.torso_top(2)], 'k-', 'LineWidth', 3);
    plot(ax, [b.hip(1) b.stance_knee(1) b.stance_foot(1)], ...
             [b.hip(2) b.stance_knee(2) b.stance_foot(2)], '-o', ...
             'Color', [0.85 0.2 0.2], 'LineWidth', 2.5, 'MarkerSize', 5);
    plot(ax, [b.hip(1) b.swing_knee(1) b.swing_foot(1)], ...
             [b.hip(2) b.swing_knee(2) b.swing_foot(2)], '--s', ...
             'Color', [0.2 0.4 0.85], 'LineWidth', 2, 'MarkerSize', 5);
    plot(ax, b.hip(1), b.hip(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    if norm(F) > 5
        quiver(ax, sh(1), sh(2), A_SCALE * F(1), A_SCALE * F(2), 0, ...
               'Color', [0.1 0.6 0.2], 'LineWidth', 2.5, 'MaxHeadSize', 0.8);
    end
    plot(ax, sh(1), sh(2), 'o', 'Color', [0.1 0.6 0.2], 'MarkerFaceColor', [0.1 0.6 0.2]);

    title(ax, {ttl, sprintf('t = %.2f s    F_x = %6.1f N    F_z = %6.1f N', ...
                             t_fr(k), F(1), F(2))});
    cur.Value = t_fr(k);
    drawnow;

    if ~isempty(gifpath)
        [im, map] = rgb2ind(frame2im(getframe(fig)), 256);
        if first
            imwrite(im, map, gifpath, 'gif', 'LoopCount', Inf, 'DelayTime', DT / speed);
            first = false;
        else
            imwrite(im, map, gifpath, 'gif', 'WriteMode', 'append', 'DelayTime', DT / speed);
        end
    else
        pause(DT / speed);
    end
end
if ~isempty(gifpath), fprintf('ch3_ub_animate: wrote %s\n', gifpath); end
end

% ---------------------------------------------------------------------------
function [z, p] = load_gait(f)
S = load(f);
pairs = {'z', 'o'; 'zC3', 'oC3'; 'zC2', 'oC2'; 'zC', 'oC'};
for i = 1:size(pairs, 1)
    if isfield(S, pairs{i, 1}) && isfield(S, pairs{i, 2})
        z = S.(pairs{i, 1});  p = S.(pairs{i, 2}).p;  return;
    end
end
error('ch3_ub_animate:file', ...
      '%s holds no {z, o} pair (fields: %s).', f, strjoin(fieldnames(S).', ', '));
end
