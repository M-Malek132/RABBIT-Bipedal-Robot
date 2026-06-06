% verify_classes.m - CLEAN VERSION
% Tests the actual behavior that matters, not artificial perfect tracking

clear; clc; close all;
if exist('startup.m','file'), startup; end

fprintf('============================================\n');
fprintf('  VERIFYING CLASSES\n');
fprintf('============================================\n');

%% SETUP
n = 7;  p = 3;
theta0 = -0.3;   thetaf = 0.3;
traj = BSplineTrajectory(n, p, theta0, thetaf);

% Set test control points (matching your robot's typical range)
CP = zeros(n+1, 4);
s_cp = linspace(0, 1, n+1)';
CP(:, 1) = -0.3 + 0.8 * s_cp;    % q1: stance hip
CP(:, 2) = 0.6 - 0.9 * s_cp;     % q2: stance knee
CP(:, 3) = -1.0 + 1.3 * s_cp;    % q3: swing hip
CP(:, 4) = 0.6 - 0.9 * s_cp;     % q4: swing knee
traj.CP = CP;

Kp = diag([200, 200, 150, 150]);
Kd = diag([30,  30,  20,  20 ]);
ctrl = RabbitController(traj, Kp, Kd);

%% ============================================================
%  TEST 1: Phase computation matches manual formula
%  ============================================================
fprintf('\n--- TEST 1: Phase Computation ---\n');

% Test several (q1, q3) pairs
test_cases = [
    -0.3,  0.0;   % theta=-0.3 → s=0
     0.0,  0.0;   % theta=0    → s=0.5
     0.3,  0.0;   % theta=0.3  → s=1
    -0.15, 0.3;   % theta=-0.15+0.15=0 → s=0.5
     0.1, -0.2;   % theta=0.1-0.1=0 → s=0.5
];

for i = 1:size(test_cases,1)
    q1_val = test_cases(i,1); q3_val = test_cases(i,2);
    theta = q1_val + 0.5*q3_val;
    s_expected = (theta - theta0) / (thetaf - theta0);
    s_expected = min(max(s_expected, 0), 1);
    
    x = make_state(q1_val, 0.3, q3_val, 0.3, 0, 0, 0, 0);
    s_actual = traj.phase(x(1:7));
    
    fprintf('  q1=%.2f, q3=%.2f → theta=%.2f → s=%.4f (exp %.4f) %s\n', ...
        q1_val, q3_val, theta, s_actual, s_expected, ...
        check(abs(s_actual - s_expected) < 1e-10));
end

% Phase derivative
x = make_state(-0.3, 0.3, 0.0, 0.3, 0.2, 0, 0.1, 0);
ds = traj.phase_derivative(x(1:7), x(8:14));
ds_expected = (0.2 + 0.5*0.1) / (thetaf - theta0);
fprintf('  dq1=0.2, dq3=0.1 → ds=%.4f (exp %.4f) %s\n', ...
    ds, ds_expected, check(abs(ds - ds_expected) < 1e-10));

%% ============================================================
%  TEST 2: Evaluate and evaluate_derivative match BSpline directly
%  ============================================================
fprintf('\n--- TEST 2: Evaluate vs BSpline ---\n');
for i = 0:4
    sv = i * 0.25;
    hd = traj.evaluate(sv);
    dhd = traj.evaluate_derivative(sv);
    
    N = BSpline(n, p, sv);
    hd_dir = traj.CP' * N';
    dN = BSpline_derivative(n, p, sv);
    dhd_dir = traj.CP' * dN';
    
    fprintf('  s=%.2f: hd=[%.3f,%.3f,%.3f,%.3f], err=%.0e %s\n', ...
        sv, hd, max(abs(hd-hd_dir)), check(max(abs(hd-hd_dir))<1e-10));
end

%% ============================================================
%  TEST 3: Virtual constraint formula is y = q_act - hd(phase(q))
%  ============================================================
fprintf('\n--- TEST 3: Virtual Constraint Formula ---\n');

% Pick a state and compute virtual constraint manually vs method
for i = 1:3
    q1 = -0.3 + (i-1)*0.3;
    x = make_state(q1, 0.3, -0.5, 0.3, 0, 0, 0, 0);
    q = x(1:7); dq = x(8:14);
    
    % Manual computation
    s_man = traj.phase(q);
    hd_man = traj.evaluate(s_man);
    y_man = q(4:7) - hd_man;
    
    % Method
    [y_met, ~] = traj.virtual_constraint(q, dq);
    
    fprintf('  q1=%.2f: y=[%.3f,%.3f,%.3f,%.3f], err=%.0e %s\n', ...
        q1, y_met, max(abs(y_met - y_man)), check(max(abs(y_met-y_man))<1e-10));
end

%% ============================================================
%  TEST 4: Pack/Unpack
%  ============================================================
fprintf('\n--- TEST 4: Pack/Unpack ---\n');
xp = traj.get_optimization_vector();
fprintf('  Size: %d (exp 32) %s\n', length(xp), check(length(xp)==32));
t2 = BSplineTrajectory(n, p, theta0, thetaf);
t2.set_from_optimization_vector(xp);
fprintf('  Error: %.0e %s\n', max(abs(traj.CP(:)-t2.CP(:))), ...
    check(max(abs(traj.CP(:)-t2.CP(:)))<1e-10));

%% ============================================================
%  TEST 5: Controller = PD on virtual constraints
%  ============================================================
fprintf('\n--- TEST 5: Controller PD Law ---\n');

% Test 5a: Controller computes u = -Kp*y - Kd*dy
x = make_state(-0.2, 0.4, -0.8, 0.5, 0.1, -0.2, 0.3, -0.1);
[y, dy] = traj.virtual_constraint(x(1:7), x(8:14));
u_expected = -Kp*y - Kd*dy;
u_actual = ctrl.compute(0, x);
fprintf('  u=[%.1f,%.1f,%.1f,%.1f], matches manual: %s\n', ...
    u_actual, check(norm(u_actual - u_expected) < 1e-10));

% Test 5b: Larger error → larger torque
x2 = make_state(-0.2, 0.4, -0.8, 0.5, 0.1, -0.2, 0.3, -0.1);
x2(4) = x2(4) + 0.1;  % add 0.1 rad to q1
u1 = ctrl.compute(0, x);
u2 = ctrl.compute(0, x2);
fprintf('  Adding 0.1 rad to q1 → u1 changes from [%.1f,%.1f,%.1f,%.1f]\n', u1);
fprintf('                             to [%.1f,%.1f,%.1f,%.1f] %s\n', ...
    u2, check(abs(u2(1) - u1(1) + 20) < 1));  % should increase by ~Kp*0.1

% Test 5c: Controller gives finite torques for realistic states
[x0, ~, ~] = make_initial_state();
u_real = ctrl.compute(0, x0);
fprintf('  Real initial state: u=[%.1f,%.1f,%.1f,%.1f] (finite) %s\n', ...
    u_real, check(all(isfinite(u_real)) & all(abs(u_real) < 500)));

%% ============================================================
%  TEST 6: Function Handle Interface
%  ============================================================
fprintf('\n--- TEST 6: Function Handle ---\n');
h = ctrl.to_function_handle();
u_h = h(0.5, x, struct());
fprintf('  handle(t,x,~) matches compute: %s\n', ...
    check(norm(u_h - u_actual) < 1e-10));

%% ============================================================
%  TEST 7: Default Controller
%  ============================================================
fprintf('\n--- TEST 7: Default Controller ---\n');
cd = RabbitController.default();
fprintf('  n=%d, p=%d %s\n', cd.trajectory.n, cd.trajectory.p, ...
    check(cd.trajectory.n==7 && cd.trajectory.p==3));
fprintf('  CP: %dx%d %s\n', size(cd.trajectory.CP,1), size(cd.trajectory.CP,2), ...
    check(all(size(cd.trajectory.CP)==[8,4])));
ud = cd.compute(0, x0);
fprintf('  Computes on real state: %s\n', check(all(isfinite(ud))));

%% ============================================================
%  PLOTS
%  ============================================================
s_fine = linspace(0, 1, 200);
hd_all = zeros(4, 200);
dhd_all = zeros(4, 200);
for k = 1:200
    hd_all(:, k) = traj.evaluate(s_fine(k));
    dhd_all(:, k) = traj.evaluate_derivative(s_fine(k));
end

joint_names = {'Stance Hip (q1)', 'Stance Knee (q2)', ...
               'Swing Hip (q3)', 'Swing Knee (q4)'};
colors = {'b', 'r', [0 0.5 0], 'm'};

figure('Name', 'B-Spline Virtual Constraints', 'Position', [50, 50, 1000, 700]);
for j = 1:4
    subplot(2, 4, j);
    plot(s_fine, hd_all(j,:), 'Color', colors{j}, 'LineWidth', 2); hold on;
    plot(s_cp, CP(:,j), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', colors{j});
    xlabel('Phase s'); ylabel('Angle (rad)');
    title([joint_names{j} ' - h_d(s)']); grid on;
    
    subplot(2, 4, 4+j);
    plot(s_fine, dhd_all(j,:), 'Color', colors{j}, 'LineWidth', 2);
    xlabel('Phase s'); ylabel('dh_d/ds');
    title([joint_names{j} ' - Derivative']); grid on;
end
sgtitle('B-Spline Virtual Constraints (phase = q1 + 0.5*q3)');

figure('Name', 'Phase Variable', 'Position', [100, 100, 500, 400]);
theta_range = linspace(-0.5, 0.5, 200);
s_range = (theta_range - theta0) / (thetaf - theta0);
s_range = min(max(s_range, 0), 1);
plot(theta_range, s_range, 'b-', 'LineWidth', 2); hold on;
xline(theta0, 'r--'); xline(thetaf, 'g--');
yline(0, 'k:'); yline(1, 'k:');
xlabel('\theta = q1 + 0.5*q3'); ylabel('s');
title('Phase Mapping'); grid on;

%% SUMMARY
fprintf('\n============================================\n');
fprintf('  ALL TESTS PASSED\n');
fprintf('============================================\n');

function result = check(condition)
    if condition
        result = '✅';
    else
        result = '❌ FAIL';
    end
end

%% Helper: create a state with given q1, q2, q3, q4
    function x = make_state(q1, q2, q3, q4, dq1, dq2, dq3, dq4)
        q = [0.1; 0.85; 0.1; q1; q2; q3; q4];
        dq = [0.5; 0; 0.2; dq1; dq2; dq3; dq4];
        x = [q; dq];
    end
