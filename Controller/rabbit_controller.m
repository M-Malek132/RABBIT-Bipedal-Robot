function u = rabbit_controller(x)
%RABBIT_CONTROLLER  Computed-torque virtual-constraint law  (ILLUSTRATIVE ONLY).
%
%   u = rabbit_controller(x)
%
% ==> NON-RUNNABLE REFERENCE. This predates the HZD pipeline and still calls
%     desired_gait / desired_gait_velocity, which DO NOT EXIST in the repo. It
%     is kept only to illustrate the computed-torque form of the virtual-
%     constraint tracking law. For the working code use, instead:
%       * hzd_control_and_dynamics.m  -- the fixed-gain PD virtual-constraint
%         law actually run by the simulator/optimizer, and
%       * Controller/rabbit_clf_qp_controller.m -- the CLF-QP alternative.
%
% Fixed here (was a silent index-convention bug): the phase variable now comes
% from the canonical theta_of_q / dtheta_dq_of (theta = qt + q1 + 0.5*q2 =
% q(3)+q(4)+0.5*q(5)) instead of the stale, torso-omitting "q(4)+0.5*q(6)".
% Keeping a second hand-written copy of the phase variable is exactly the kind
% of drift that breaks the leg-relabeling conventions -- always call the shared
% helpers.

    % 1. Extract States
    nq = 7;
    q  = x(1:nq);
    dq = x(nq+1:end);

    % 2. Get Robot Dynamics Matrices
    % You need these from your model (e.g., via Symbolic Toolbox or Euler-Lagrange)
    M_matrix = M(q);       % 7x7 Mass Matrix
    C = V([q; dq]); % 7x1 Coriolis/Gravity vector

    % 3. Define Phase Variable (s) -- use the shared, canonical helpers so this
    % cannot drift from the rest of the pipeline (theta = c*q, c from
    % dtheta_dq_of). theta0/thetaf are the illustrative phase endpoints
    % (thetaf mirrors p.theta_plus = 0.3).
    theta  = theta_of_q(q);
    theta0 = -0.3; thetaf = 0.3;
    ds_dtheta = 1/(thetaf - theta0);
    s = (theta - theta0) * ds_dtheta;
    s = min(max(s,0),1);
    dtheta_dt = dtheta_dq_of(q) * dq;          % theta_dot = c*dq

    % 4. Virtual Constraints
    % hd should be the desired joint positions; d(hd)/dt via the chain rule
    % d(hd)/ds * ds/dtheta * dtheta/dt (the ds/dtheta factor was missing before).
    hd  = desired_gait(s);
    dhd = desired_gait_velocity(s) * ds_dtheta * dtheta_dt;

    % Tracking Error
    y  = q(4:7) - hd;
    dy = dq(4:7) - dhd;

    % 5. Computed Torque Control
    % We want ddq_desired = -Kp*y - Kd*dy
    % So, u = M * (-Kp*y - Kd*dy) + C
    Kp = diag([400, 400, 300, 300]);
    Kd = diag([40, 40, 30, 30]);

    % We assume only joints 4,5,6,7 are actuated
    % Map the control input to the actuated joints
    B = [zeros(3,4); eye(4)]; 
    
    % Compute required acceleration for tracking
    ddq_ref = -Kp*y - Kd*dy;
    
    % The control law
    u = M_matrix(4:7, :) * [zeros(3,1); ddq_ref] + C(4:7);
end
