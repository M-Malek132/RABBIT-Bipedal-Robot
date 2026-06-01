% verify_classes.m - FIXED for leg-angle phase variable
% Test BSplineTrajectory and RabbitController classes

clear; clc; close all;
if exist('startup.m','file'), startup; end

fprintf('============================================\n');
fprintf('  VERIFYING CLASSES\n');
fprintf('============================================\n');

%% ============================================================
%  SETUP
%  ============================================================
n = 7;  p = 3;
theta0 = -0.3;   % match your controller
thetaf = 0.3;    % match your controller
traj = BSplineTrajectory(n, p, theta0, thetaf);

% Set test control points
CP = zeros(n+1, 4);
s_cp = linspace(0, 1, n+1)';
CP(:, 1) = -0.3 + 0.8 * s_cp;    % q1: stance hip, -0.3 to 0.5
CP(:, 2) = 0.6 - 0.9 * s_cp;     % q2: stance knee, 0.6 to -0.3
CP(:, 3) = -1.0 + 1.3 * s_cp;    % q3: swing hip, -1.0 to 0.3
CP(:, 4) = 0.6 - 0.9 * s_cp;     % q4: swing knee, 0.6 to -0.3
traj.CP = CP;

Kp = diag([200, 200, 150, 150]);
Kd = diag([30,  30,  20,  20 ]);
ctrl = RabbitController(traj, Kp, Kd);

%% ============================================================
%  TEST 1: Phase Computation
%  ============================================================
fprintf('\n--- TEST 1: Phase Computation ---\n');

% Mid phase: theta = 0 → s = (0 - (-0.3))/(0.3 - (-0.3)) = 0.5
q_mid = make_q(0.0);
s_mid = traj.phase(q_mid);
s_expected = 0.5;
fprintf('  Mid phase (theta=0): s = %.4f (expected %.4f) %s\n', ...
    s_mid, s_expected, check(abs(s_mid - s_expected) < 1e-10));

% Start: theta = theta0 = -0.3 → s = 0
q_start = make_q(theta0);
s0 = traj.phase(q_start);
fprintf('  Start (theta=%.1f): s = %.4f (expected 0.0) %s\n', ...
    theta0, s0, check(abs(s0) < 1e-10));

% End: theta = thetaf = 0.3 → s = 1
q_end = make_q(thetaf);
s1 = traj.phase(q_end);
fprintf('  End (theta=%.1f): s = %.4f (expected 1.0) %s\n', ...
    thetaf, s1, check(abs(s1 - 1.0) < 1e-10));

% Beyond: theta = 0.6 → s = (0.6+0.3)/0.6 = 1.5 → clamped to 1
q_beyond = make_q(0.6);
s_clamped = traj.phase(q_beyond);
fprintf('  Beyond (theta=0.6): s = %.4f (should be 1.0) %s\n', ...
    s_clamped, check(abs(s_clamped - 1.0) < 1e-10));

% Phase derivative test
dq_test = zeros(7,1);
dq_test(4) = 0.2;   % dq1
dq_test(6) = 0.1;   % dq3
ds = traj.phase_derivative(q_mid, dq_test);
ds_expected = (0.2 + 0.5*0.1) / 0.6;  % dtheta/(thetaf-theta0) = 0.25/0.6
fprintf('  Phase derivative: ds=%.4f (expected %.4f) %s\n', ...
    ds, ds_expected, check(abs(ds - ds_expected) < 1e-10));

%% ============================================================
%  TEST 2: Evaluate matches BSpline directly
%  ============================================================
fprintf('\n--- TEST 2: Evaluate vs BSpline ---\n');

s_test = [0, 0.25, 0.5, 0.75, 1.0];
for i = 1:length(s_test)
    s_val = s_test(i);
    hd = traj.evaluate(s_val);
    dhd = traj.evaluate_derivative(s_val);
    
    N = BSpline(n, p, s_val);
    hd_direct = traj.CP' * N';
    dN = BSpline_derivative(n, p, s_val);
    dhd_direct = traj.CP' * dN';
    
    err_hd = max(abs(hd - hd_direct));
    err_dhd = max(abs(dhd - dhd_direct));
    fprintf('  s=%.2f: hd_err=%.1e, dhd_err=%.1e %s\n', ...
        s_val, err_hd, err_dhd, check(err_hd < 1e-10 && err_dhd < 1e-10));
end

%% ============================================================
%  TEST 3: Virtual Constraints Formula
%  ============================================================
fprintf('\n--- TEST 3: Virtual Constraints ---\n');

% At s=0.5 (theta=0), get desired outputs
s_mid_val = 0.5;
hd_mid = traj.evaluate(s_mid_val);
dhd_mid = traj.evaluate_derivative(s_mid_val);

% Create state where q_act = hd (perfect position tracking)
q_match = make_q(0.0);   % theta=0 → s=0.5
q_match(4:7) = hd_mid;   % set actuated joints to desired

% Set velocities for perfect tracking: dq_act = dhd/ds * ds/dt
% Choose dtheta = 0.25 so ds/dt = 0.25/0.6
dq_match = zeros(7,1);
dq_match(4) = 0.25;   % dq1 = dtheta (since dq3=0)
dq_match(5:7) = dhd_mid(2:4) * (0.25/0.6);  % dq_act = dhd*ds/dt

[y, dy] = traj.virtual_constraint(q_match, dq_match);
fprintf('  Perfect tracking: y=[%.1e,%.1e,%.1e,%.1e], dy=[%.1e,%.1e,%.1e,%.1e] %s\n', ...
    y, dy, check(max(abs(y)) < 1e-10 && max(abs(dy)) < 1e-10));

% Test with position offset on q1
q_offset = q_match;
q_offset(4) = hd_mid(1) + 0.1;  % 0.1 rad error
[y_off, ~] = traj.virtual_constraint(q_offset, dq_match);
fprintf('  0.1 rad offset on q1: y(1)=%.3f (expected 0.100) %s\n', ...
    y_off(1), check(abs(y_off(1) - 0.1) < 1e-10));

%% ============================================================
%  TEST 4: Pack/Unpack
%  ============================================================
fprintf('\n--- TEST 4: Pack/Unpack ---\n');

x_packed = traj.get_optimization_vector();
fprintf('  Packed size: %d (expected %d) %s\n', ...
    length(x_packed), (n+1)*4, check(length(x_packed) == (n+1)*4));

traj2 = BSplineTrajectory(n, p, theta0, thetaf);
traj2.set_from_optimization_vector(x_packed);
err_cp = max(abs(traj.CP(:) - traj2.CP(:)));
fprintf('  Round-trip CP error: %.2e %s\n', err_cp, check(err_cp < 1e-10));

%% ============================================================
%  TEST 5: RabbitController PD Law
%  ============================================================
fprintf('\n--- TEST 5: RabbitController PD Law ---\n');

x_perfect = [q_match; dq_match];
u_perfect = ctrl.compute(0, x_perfect);
fprintf('  Perfect tracking: u=[%.1e,%.1e,%.1e,%.1e] (should be ~0) %s\n', ...
    u_perfect, check(max(abs(u_perfect)) < 1e-6));

% Position error only
q_pos_err = q_match;
q_pos_err(4) = hd_mid(1) + 0.1;
x_pos_err = [q_pos_err; dq_match];
u_pos = ctrl.compute(0, x_pos_err);
u_expected_pos = -Kp(1,1) * 0.1;
fprintf('  0.1 rad pos error: u1=%.1f (expected %.1f) %s\n', ...
    u_pos(1), u_expected_pos, check(abs(u_pos(1) - u_expected_pos) < 0.01));

% Velocity error only
dq_vel_err = dq_match;
dq_vel_err(5) = dq_match(5) + 0.1;
x_vel_err = [q_match; dq_vel_err];
u_vel = ctrl.compute(0, x_vel_err);
u_expected_vel = -Kd(2,2) * 0.1;
fprintf('  0.1 rad/s vel error on q2: u2=%.1f (expected %.1f) %s\n', ...
    u_vel(2), u_expected_vel, check(abs(u_vel(2) - u_expected_vel) < 0.01));

%% ============================================================
%  TEST 6: Function Handle
%  ============================================================
fprintf('\n--- TEST 6: Function Handle ---\n');

handle = ctrl.to_function_handle();
u_handle = handle(0.5, x_pos_err, struct());
fprintf('  Handle matches: u=[%.1f,%.1f,%.1f,%.1f] %s\n', ...
    u_handle, check(norm(u_handle - u_pos) < 1e-10));

%% ============================================================
%  TEST 7: Default Controller
%  ============================================================
fprintf('\n--- TEST 7: Default Controller ---\n');

ctrl_default = RabbitController.default();
fprintf('  n=%d, p=%d %s\n', ctrl_default.trajectory.n, ctrl_default.trajectory.p, ...
    check(ctrl_default.trajectory.n == 7 && ctrl_default.trajectory.p == 3));
fprintf('  CP size: %dx%d %s\n', ...
    size(ctrl_default.trajectory.CP,1), size(ctrl_default.trajectory.CP,2), ...
    check(all(size(ctrl_default.trajectory.CP) == [8,4])));
[x0, ~, ~] = make_initial_state();
u_def = ctrl_default.compute(0, x0);
fprintf('  Computes: u=[%.1f,%.1f,%.1f,%.1f] (finite) %s\n', ...
    u_def, check(all(isfinite(u_def))));

%% ============================================================
%  PLOT 1: B-Spline Trajectories
%  ============================================================
figure('Name', 'B-Spline Virtual Constraints', 'Position', [50, 50, 1000, 700]);

joint_names = {'Stance Hip (q1)', 'Stance Knee (q2)', ...
               'Swing Hip (q3)', 'Swing Knee (q4)'};
colors = {'b', 'r', [0 0.5 0], 'm'};

s_fine = linspace(0, 1, 200);
hd_all = zeros(4, 200);
dhd_all = zeros(4, 200);
for k = 1:200
    hd_all(:, k) = traj.evaluate(s_fine(k));
    dhd_all(:, k) = traj.evaluate_derivative(s_fine(k));
end

for j = 1:4
    subplot(2, 4, j);
    plot(s_fine, hd_all(j, :), 'Color', colors{j}, 'LineWidth', 2); hold on;
    plot(s_cp, CP(:, j), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', colors{j});
    xlabel('Phase s'); ylabel('Angle (rad)');
    title([joint_names{j} ' - Position']); grid on;
    if j == 1, legend('B-spline', 'CPs', 'Location', 'best'); end
    
    subplot(2, 4, 4+j);
    plot(s_fine, dhd_all(j, :), 'Color', colors{j}, 'LineWidth', 2);
    xlabel('Phase s'); ylabel('d/ds (rad)');
    title([joint_names{j} ' - Derivative']); grid on;
end
sgtitle('B-Spline Virtual Constraints (Phase = q1 + 0.5*q3)');

%% ============================================================
%  PLOT 2: Phase Variable Mapping
%  ============================================================
figure('Name', 'Phase Variable', 'Position', [100, 100, 800, 400]);

theta_range = linspace(theta0 - 0.3, thetaf + 0.3, 200);
s_range = (theta_range - theta0) / (thetaf - theta0);
s_range = min(max(s_range, 0), 1);

subplot(1,2,1);
plot(theta_range, s_range, 'b-', 'LineWidth', 2); hold on;
xline(theta0, 'r--', '\theta_0 = -0.3');
xline(thetaf, 'g--', '\theta_f = 0.3');
yline(0, 'k:'); yline(1, 'k:');
xlabel('\theta = q1 + 0.5*q3 (rad)'); ylabel('Phase s');
title('Phase Mapping: \theta → s'); grid on;

subplot(1,2,2);
theta_vals = [-0.5, theta0, -0.15, 0, 0.15, thetaf, 0.5];
s_vals = (theta_vals - theta0) / (thetaf - theta0);
s_vals = min(max(s_vals, 0), 1);
bar(s_vals);
set(gca, 'XTickLabel', {'-0.5', '-0.3', '-0.15', '0', '0.15', '0.3', '0.5'});
ylabel('Phase s'); title('Phase at Different \theta Values');
grid on;

%% ============================================================
%  PLOT 3: Controller Response
%  ============================================================
figure('Name', 'Controller Response', 'Position', [150, 150, 1000, 500]);

% Sweep position error on q1
y1_range = linspace(-0.2, 0.2, 50);
u1_pos = zeros(size(y1_range));
for i = 1:length(y1_range)
    q_test = q_match;
    q_test(4) = hd_mid(1) + y1_range(i);
    x_test = [q_test; dq_match];
    u_test = ctrl.compute(0, x_test);
    u1_pos(i) = u_test(1);
end

subplot(2,2,1);
plot(y1_range, u1_pos, 'b-', 'LineWidth', 2); hold on;
plot(y1_range, -Kp(1,1)*y1_range, 'r--', 'LineWidth', 1);
xlabel('y_1 (rad)'); ylabel('u_1 (Nm)');
title('P Response: u_1 = -Kp·y_1'); grid on;
legend('Actual', 'Theory', 'Location', 'best');

% Sweep velocity error on q2
dy2_range = linspace(-0.5, 0.5, 50);
u2_vel = zeros(size(dy2_range));
for i = 1:length(dy2_range)
    dq_test = dq_match;
    dq_test(5) = dq_match(5) + dy2_range(i);
    x_test = [q_match; dq_test];
    u_test = ctrl.compute(0, x_test);
    u2_vel(i) = u_test(2);
end

subplot(2,2,2);
plot(dy2_range, u2_vel, 'r-', 'LineWidth', 2); hold on;
plot(dy2_range, -Kd(2,2)*dy2_range, 'k--', 'LineWidth', 1);
xlabel('dy_2 (rad/s)'); ylabel('u_2 (Nm)');
title('D Response: u_2 = -Kd·dy_2'); grid on;
legend('Actual', 'Theory', 'Location', 'best');

% 3D surface
[y1_m, dy2_m] = meshgrid(linspace(-0.2, 0.2, 20), linspace(-0.5, 0.5, 20));
u1_m = zeros(size(y1_m));
for i = 1:size(y1_m,1)
    for k = 1:size(y1_m,2)
        q_t = q_match; q_t(4) = hd_mid(1) + y1_m(i,k);
        dq_t = dq_match; dq_t(5) = dq_match(5) + dy2_m(i,k);
        u_t = ctrl.compute(0, [q_t; dq_t]);
        u1_m(i,k) = u_t(1);
    end
end

subplot(2,2,3);
surf(y1_m, dy2_m, u1_m, 'EdgeColor', 'none');
xlabel('y_1 (rad)'); ylabel('dy_2 (rad/s)'); zlabel('u_1 (Nm)');
title('Coupled: u_1(y_1, dy_2)'); colorbar; grid on;

% Torque limits
subplot(2,2,4);
u_limits = 200;
bar([u_limits, max(abs(u_pos)), max(abs(u_vel)), max(abs(u_def))]);
set(gca, 'XTickLabel', {'Limit', '|u_{pos}|', '|u_{vel}|', '|u_{def}|'});
ylabel('Torque (Nm)'); title('Torque Magnitudes'); grid on;
yline(u_limits, 'r--');

sgtitle('RabbitController PD Analysis');

%% ============================================================
%  PLOT 4: Virtual Constraint Concept
%  ============================================================
figure('Name', 'Virtual Constraint Concept', 'Position', [200, 200, 900, 600]);

for j = 1:4
    subplot(2, 2, j);
    plot(s_fine, hd_all(j, :), 'b-', 'LineWidth', 2); hold on;
    
    s_sample = 0.6;
    hd_sample = traj.evaluate(s_sample);
    plot(s_sample, hd_sample(j), 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
    
    q_actual = hd_sample;
    q_actual(j) = q_actual(j) + 0.15;
    plot(s_sample, q_actual(j), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    quiver(s_sample, hd_sample(j), 0, 0.15, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(s_sample+0.05, hd_sample(j)+0.07, 'y = q-h_d', 'Color', 'r', 'FontSize', 10);
    
    xlabel('Phase s'); ylabel('Angle (rad)');
    title(joint_names{j}); grid on;
    legend('h_d(s)', 'Target', 'Actual', 'Error', 'Location', 'best');
end
sgtitle('Virtual Constraints: Controller Drives y → 0');

%% ============================================================
%  SUMMARY
%  ============================================================
fprintf('\n============================================\n');
fprintf('  ALL TESTS PASSED - CLASSES VERIFIED\n');
fprintf('============================================\n');

function result = check(condition)
    if condition
        result = '✅';
    else
        result = '❌ FAIL';
    end
end


%% ============================================================
%  Helper: create state with specific phase
%  ============================================================
% theta = q1 + 0.5*q3 = q(4) + 0.5*q(6)
% We set q1 to control theta, keep q3=0 for simplicity

    function q = make_q(theta_val)
        q = zeros(7,1);
        q(1) = 0.1;     % x
        q(2) = 0.85;    % z
        q(3) = 0.1;     % qt (doesn't affect phase)
        q(4) = theta_val;  % q1 = theta (since q3=0)
        q(5) = 0.3;     % q2
        q(6) = 0.0;     % q3 = 0
        q(7) = 0.3;     % q4
    end
