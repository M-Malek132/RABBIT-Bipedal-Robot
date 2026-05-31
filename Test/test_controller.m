% quick_controller_check.m

clear; clc;

%% Initialize
params = init_bspline_params();

%% Create a test state
% Mid-step configuration
q_test = [0.2;      % x: moved forward a bit
          0.85;     % z: hip height
          0.1;      % qt: torso tilted forward
          -0.1;     % q1: stance knee extending
          0.2;      % q2: stance hip
          0.0;      % q3: swing knee
          0.2];     % q4: swing hip

dq_test = [1.0;     % dx: walking forward
           0.0;     % dz
           0.25;    % dqt
           -0.8;    % dq1
           1.2;     % dq2
           -0.6;    % dq3
           0.9];    % dq4

x_test = [q_test; dq_test];

%% Call controller
u = rabbit_controller_bspline(0.5, x_test, params);

%% Display results
fprintf('Controller Test:\n');
fprintf('  Phase (qt) = %.3f (range: [%.2f, %.2f])\n', q_test(3), params.theta0, params.thetaf);
fprintf('  s = %.3f\n', (q_test(3) - params.theta0)/(params.thetaf - params.theta0));
fprintf('  Torques: [%.2f, %.2f, %.2f, %.2f] Nm\n', u);

%% Check virtual constraint error
s = (q_test(3) - params.theta0) / (params.thetaf - params.theta0);
s = min(max(s, 0), 1);
hd = desired_gait_bspline(s, params.ControlPoints, params.n, params.p);
y = q_test(4:7) - hd;
fprintf('  Virtual constraint errors: [%.3f, %.3f, %.3f, %.3f]\n', y);