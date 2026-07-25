function fig = ch3_plot_gait(z, p, savepath)
%CH3_PLOT_GAIT  Six-panel summary of an optimized gait.
%
%   fig = ch3_plot_gait(z, p)
%   fig = ch3_plot_gait(z, p, 'Results/ch3_gait.png')
%
% Panels:
%   1  stick-figure frames through the step (the sanity check that beats any
%      number: does it look like walking?)
%   2  virtual constraints -- actual joint angles vs the commanded Bezier
%      profiles, against the phase s
%   3  joint torques
%   4  ground reaction force, with the friction cone and Fz limits drawn
%   5  the transverse variables eta = [y; ydot] across the step
%   6  swing-foot height, with the clearance requirement drawn
%
% Requires a normal (JVM-enabled) MATLAB session; solves can run under -nojvm
% but plotting cannot.
%
% See also CH3_REPORT, CH3_ANIMATE.

E = ch3_col_eval(z, p);
[X, T, alpha] = ch3_col_unpack(z, p);
N = E.N;
t = linspace(0, T, N);

fig = figure('Color','w','Position',[80 80 1250 760]);

%% 1. stick figure
subplot(2,3,1); hold on; grid on; axis equal
idx = round(linspace(1, N, 6));
cols = parula(numel(idx));
for j = 1:numel(idx)
    q = X(1:p.nq, idx(j));
    b = ch3_body_points(q);
    c = cols(j,:);
    plot([b.hip(1) b.torso_top(1)],   [b.hip(2) b.torso_top(2)],   '-', 'Color', c, 'LineWidth', 2);
    plot([b.hip(1) b.stance_knee(1) b.stance_foot(1)], ...
         [b.hip(2) b.stance_knee(2) b.stance_foot(2)], '-o', 'Color', c, 'MarkerSize', 3);
    plot([b.hip(1) b.swing_knee(1) b.swing_foot(1)], ...
         [b.hip(2) b.swing_knee(2) b.swing_foot(2)], '--s', 'Color', c, 'MarkerSize', 3);
end
yline(0, 'k-', 'LineWidth', 1.5);
title(sprintf('gait: %.3f m in %.3f s (%.3f m/s)', E.L_step, T, E.L_step/T));
xlabel('x [m]'); ylabel('z [m]');

%% 2. virtual constraints
subplot(2,3,2); hold on; grid on
s_nodes = zeros(1,N); Y = zeros(p.ny, N); YD = zeros(p.ny, N);
for k = 1:N
    s_nodes(k) = ch3_phase(X(:,k), p);
    Y(:,k)  = p.H * X(1:p.nq,k);
    YD(:,k) = ch3_yd(alpha, s_nodes(k), p);
end
co = lines(p.ny);
for i = 1:p.ny
    plot(s_nodes, Y(i,:),  '-',  'Color', co(i,:), 'LineWidth', 1.6);
    plot(s_nodes, YD(i,:), '--', 'Color', co(i,:), 'LineWidth', 1.0);
end
xlabel('phase s'); ylabel('joint angle [rad]');
title('virtual constraints: q (solid) vs y_d(s) (dashed)');

%% 3. torques
subplot(2,3,3); hold on; grid on
plot(t, E.u.', 'LineWidth', 1.4);
if p.limits.enable.torque
    yline( p.limits.u_max, 'r--'); yline(-p.limits.u_max, 'r--');
end
xlabel('t [s]'); ylabel('u [Nm]');
title(sprintf('joint torques (peak %.1f Nm)', max(abs(E.u(:)))));
legend({'st hip','st knee','sw hip','sw knee'}, 'Location','best','FontSize',7);

%% 4. ground reaction force
subplot(2,3,4); hold on; grid on
plot(t, E.lam(1,:), 'LineWidth', 1.4);
plot(t, E.lam(2,:), 'LineWidth', 1.4);
plot(t, p.limits.mu_s * E.lam(2,:), 'k--');
plot(t, -p.limits.mu_s * E.lam(2,:), 'k--');
yline(p.limits.Fz_min, 'r:');
yline(0, 'k-');
xlabel('t [s]'); ylabel('force [N]');
title(sprintf('GRF (F_z min %.0f N)', min(E.lam(2,:))));
legend({'F_x','F_z','\pm\mu_s F_z'}, 'Location','best','FontSize',7);

%% 5. transverse variables
subplot(2,3,5); hold on; grid on
eta = zeros(2*p.ny, N);
for k = 1:N
    [yk, ydk] = ch3_outputs(X(:,k), alpha, p);
    eta(:,k) = [yk; ydk];
end
semilogy(t, max(abs(eta),[],1), 'LineWidth', 1.6);
set(gca,'YScale','log');
xlabel('t [s]'); ylabel('max |\eta|');
title('distance from the zero dynamics surface');

%% 6. swing-foot height
subplot(2,3,6); hold on; grid on
plot(t, E.sw_h, 'LineWidth', 1.6);
yline(p.limits.clearance, 'r--');
yline(0, 'k-');
xlabel('t [s]'); ylabel('swing-foot height [m]');
title('foot clearance');

sgtitle(sprintf('Chapter 3 gait  |  J = %.1f, %s basis, %s', ...
        ch3_col_cost(z,p), p.basis, p.controller));

if nargin >= 3 && ~isempty(savepath)
    exportgraphics(fig, savepath, 'Resolution', 150);
    fprintf('Wrote %s\n', savepath);
end

end
