function figs = ch5_plot_pendulum(runs, savedir)
%CH5_PLOT_PENDULUM  Figure 5.4: the relative-degree-4 nonlinear validation.
%
%   figs = ch5_plot_pendulum(runs)
%   figs = ch5_plot_pendulum(runs, savedir)
%
%   fig 1  Fig. 5.4 itself, one column per run:
%            row 1  control outputs y1, y2
%            row 2  control inputs tau1, tau2
%            row 3  the end effector's path in the plane, against p2min
%   fig 2  the barrier record: py2(t) against its limit, and y_rb(t)
%
% ---------------------------------------------------------- how to read fig 1
% ROW 3 IS THE FIGURE. The workspace trace shows where the end effector was
% allowed to go, and the red line is where it was not. In the baseline column
% the trace sweeps a full circle of radius 2 and dives to py = -2; in the ECBF
% columns it flattens against the red line and slides along it. That sliding is
% the barrier being ACTIVE rather than merely satisfied.
%
% ROW 1 EXPLAINS THE COST. y2 = theta2 - theta2d is what the CLF wants at zero,
% and in the ECBF columns it is driven to roughly -2 rad in the middle of the
% run -- the arm folding. The CLF is not being satisfied there and it is not
% supposed to be: the barrier row is hard, the CLF row is slacked, and this is
% what that asymmetry looks like in a trajectory. Comparing the columns of row
% 1 without row 3 next to them would read as the controller doing worse, when
% it is doing something different and harder.
%
% The torque axis in row 2 is shared across columns, for the same reason as in
% ch5_plot_springmass: cross-column comparison is the point.
%
% See also CH5_PLOT_SPRINGMASS, CH5_ANIMATE_PENDULUM, CH5_REPORT.

if nargin < 2, savedir = ''; end

n = numel(runs);
figs = gobjects(0);

%% ============================================================ fig 1: Fig 5.4
f1 = figure('Name', 'ch5: Fig 5.4 two-link pendulum, elastic actuators', ...
            'Position', [50 50 400*n 780]);

u_lim = ch5_robust_ulim(runs);

for i = 1:n
    s    = runs(i).sim;
    thd  = s.p.plant.thetad(:);
    lvl  = s.p.constraint.value;
    reach = sum(s.p.plant.l);

    % --- row 1: control outputs
    subplot(3, n, i); hold on; grid on;
    plot(s.t, s.x(1,:) - thd(1), 'LineWidth', 1.5);
    plot(s.t, s.x(2,:) - thd(2), 'LineWidth', 1.5);
    yline(0, 'k:');
    xlim([0 s.t(end)]);
    ylabel('Control outputs (rad)');
    title(runs(i).label, 'Interpreter', 'none', 'FontWeight', 'normal');
    if i == 1
        legend({'y_1 = \theta_1 - \theta_{1d}', 'y_2 = \theta_2 - \theta_{2d}'}, ...
               'Location', 'southeast');
    end

    % --- row 2: control inputs, SHARED axis
    %
    % Drawn thin because the ECBF traces CHATTER. That is not noise to be
    % filtered out of the plot: while the barrier row is active the controller
    % is holding the state on the constraint boundary with a zero-order hold,
    % the state drifts off between samples, and the next sample corrects it --
    % a discrete-time sliding mode. It is the honest picture of what enforcing
    % a hard constraint at a finite rate looks like, and the %act column in
    % ch5_report says how much of the run is spent doing it.
    ax2 = subplot(3, n, n+i); hold on; grid on;
    plot(s.t, s.u(1,:), 'LineWidth', 0.6);
    plot(s.t, s.u(2,:), 'LineWidth', 0.6);
    ylim([-u_lim u_lim]);
    xlim([0 s.t(end)]);
    ch5_note_clipping(ax2, s, u_lim, 'Nm');
    xlabel('Time (s)'); ylabel('Control inputs (Nm)');
    if i == 1, legend({'\tau_1','\tau_2'}, 'Location','southeast'); end

    % --- row 3: the end effector in the plane
    subplot(3, n, 2*n+i); hold on; grid on; axis equal;
    P   = ch5_pend_points(s.x(1:2,:), s.p);
    tip = squeeze(P(:,3,:));

    % the reachable circle, for scale
    a = linspace(0, 2*pi, 200);
    plot(reach*cos(a), reach*sin(a), ':', 'Color', [0.7 0.7 0.7]);

    plot(tip(1,:), tip(2,:), 'LineWidth', 1.3, 'Color', [0 0.45 0.74]);
    yline(lvl, 'r-', 'LineWidth', 1.5);
    plot(tip(1,1),   tip(2,1),   'ko', 'MarkerFaceColor','g', 'MarkerSize', 7);
    plot(tip(1,end), tip(2,end), 'ks', 'MarkerFaceColor','r', 'MarkerSize', 7);

    % the arm at the instant the constraint was tightest -- the pose the
    % barrier forced, which is otherwise invisible in a trace of the tip alone
    [~, k] = min(s.h);
    plot(squeeze(P(1,:,k)), squeeze(P(2,:,k)), '-o', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.6, 'MarkerSize', 4, ...
         'MarkerFaceColor', [0.2 0.2 0.2]);

    xlim([-1.1 1.1]*reach); ylim([-1.1 1.1]*reach);
    xlabel('x (m)'); ylabel('y (m)');
    if i == 1
        legend({'reach','p_2 path','p_{2min}','start','end','pose at min h'}, ...
               'Location','southoutside', 'NumColumns', 3);
    end
end

figs(end+1) = f1;

%% ====================================================== fig 2: barrier record
f2 = figure('Name', 'ch5: pendulum barrier record', 'Position', [70 70 900 520]);

subplot(2,1,1); hold on; grid on;
for i = 1:n
    s = runs(i).sim;
    plot(s.t, s.h + s.p.constraint.value, 'LineWidth', 1.5);
end
for i = 1:n
    yline(runs(i).sim.p.constraint.value, 'r--', 'LineWidth', 1.0);
end
ylabel('p^y_2   [m]');
title('End-effector height against its limits');
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
    names = {'ch5_pendulum_fig54', 'ch5_pendulum_barrier'};
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(savedir, [names{k} '.png']), ...
                       'Resolution', 150);
    end
end

end
