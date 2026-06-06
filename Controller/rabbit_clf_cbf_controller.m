function [u, delta_clf, delta_cbf] = rabbit_clf_cbf_controller(t, x, params, target_stone_idx)
% ---------------------------------------------------------------------
% rabbit_clf_cbf_controller: Robust HOCBF with CLF-CBF slack
% ---------------------------------------------------------------------
p = packParameters(params);
q = x(1:7,1);
dq = x(8:14,1);

D = D_matrix(q,p);
C = C_vector(q,dq,p);
G = G_vector(q,p);
H_dyn = C + G;
B = input_matrix();
D_inv = inv(D);

% ---------------------------------------------------------------------
% 1. CLF Formulation (Trajectory Tracking with Obstacle Bias)
% ---------------------------------------------------------------------
[y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq);
[y_d, dy_d, ddy_d] = get_ModifiedDesiredOutputs(t, q, params, target_stone_idx);

e = y - y_d;
de = J_y * dq - dy_d;
eta = [e; de];

Lfe_drift = Jdot_y_dq - J_y * D_inv * H_dyn - ddy_d;
Lge_ctrl  = J_y * D_inv * B;

F_drift = [de; Lfe_drift];
G_ctrl  = [zeros(length(e), size(B,2)); Lge_ctrl];

Kp_mat = 100 * eye(length(e));
Kd_mat = 20 * eye(length(e));
P = [Kp_mat + 0.5*Kd_mat^2, 0.5*Kd_mat;
    0.5*Kd_mat,            0.5*eye(length(e))];

V = eta' * P * eta;
LfV = 2 * eta' * P * F_drift;
LgV = 2 * eta' * P * G_ctrl;
c_clf = 10;

% ---------------------------------------------------------------------
% 2. HOCBF Formulation
% ---------------------------------------------------------------------
stone = params.stones(target_stone_idx, :);
stone_width = stone(2) - stone(1);
center_x = (stone(1) + stone(2)) / 2;

[~, p_swing, ~, ~, ~, ~] = rabbit_kinematics(q,p);
J_sw = J_swing(q,p);
Jdot_sw_dq = Jdotdq_swing(q,dq,p);

d_horizontal = p_swing(1) - center_x;
A_barrier = 0.10;
sigma = stone_width / 2;

exp_term = exp(-(d_horizontal^2) / (2 * sigma^2));
z_boundary = A_barrier * (1 - exp_term);

dz_bound_dx = A_barrier * (d_horizontal / sigma^2) * exp_term;
d2z_bound_dx2 = A_barrier * exp_term * ((1/sigma^2) - (d_horizontal^2 / sigma^4));

h = p_swing(2) - z_boundary;
Jh = J_sw(2,:) - dz_bound_dx * J_sw(1,:);
h_dot = Jh * dq;

gamma_1 = 1;
gamma_2 = 2;
psi_1 = h_dot + gamma_1 * h;

v_swing_x = J_sw(1,:) * dq;
a_drift_x = Jdot_sw_dq(1);
a_drift_z = Jdot_sw_dq(2);

h_drift = a_drift_z - (d2z_bound_dx2 * v_swing_x^2 + dz_bound_dx * a_drift_x);

Lf_psi1 = h_drift - Jh * D_inv * H_dyn + gamma_1 * h_dot;
Lg_psi1 = Jh * D_inv * B;

% ---------------------------------------------------------------------
% 3. Quadratic Program (Softened CLF + CBF Constraints)
% ---------------------------------------------------------------------
% Decision variables:
% z = [u1; u2; u3; u4; delta_clf; delta_cbf]
weight_u = 1;
weight_delta_clf = 1000;
weight_delta_cbf = 1000;

H_qp = diag([weight_u*ones(4,1); weight_delta_clf; weight_delta_cbf]);
f_qp = zeros(6,1);

% CLF:
% LfV + LgV*u + c_clf*V <= delta_clf
% => [LgV, -1, 0] * z <= -LfV - c_clf*V
A_clf = [LgV, -1, 0];
b_clf = -LfV - c_clf * V;

% HOCBF:
% Lf_psi1 + Lg_psi1*u >= -gamma_2*psi_1 - delta_cbf
% => -Lg_psi1*u - delta_cbf <= Lf_psi1 + gamma_2*psi_1
A_cbf = [-Lg_psi1, 0, -1];
b_cbf = Lf_psi1 + gamma_2 * psi_1;

A_ineq = [A_clf;
    A_cbf];
b_ineq = [b_clf;
    b_cbf];

u_max = 150;
lb = [-u_max*ones(4,1); 0; 0];
ub = [ u_max*ones(4,1); Inf; Inf];
% ... (A_ineq, b_ineq, lb, ub are defined as above)

options = optimoptions('quadprog', 'Display', 'off');
[z_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);

if exitflag == 1
    u = z_opt(1:4);
    delta_clf = z_opt(5);
    delta_cbf = z_opt(6);
else
    % --- PD FALLBACK ---
    % In case of QP infeasibility, track the modified trajectory
    % using a simple PD controller.

    warning('QP infeasible at t=%.3f. Falling back to PD control.', t);

    % Gains (Tuned to be stiff enough to maintain stability)
    Kp_pd = 150;
    Kd_pd = 20;

    % Calculate error
    error = y_d - y;
    d_error = dy_d - (J_y * dq);

    % PD Control Law
    u = Kp_pd * error + Kd_pd * d_error;

    % Apply torque saturation to match QP limits
    u = max(min(u, u_max), -u_max);

    % Slack variables set to zero as we are abandoning the constraints
    delta_clf = 0;
    delta_cbf = 0;
end

end

function [y_d_mod, dy_d_mod, ddy_d_mod] = get_ModifiedDesiredOutputs(t, q, params, target_stone_idx)
[y_d, dy_d, ddy_d] = get_DesiredOutputs(t);
p = packParameters(params);
[~, p_swing, ~, ~, ~, ~] = rabbit_kinematics(q, p);

stone = params.stones(target_stone_idx, :);
center_x = (stone(1) + stone(2)) / 2;
d_x = p_swing(1) - center_x;
sigma = (stone(2) - stone(1)) / 2;

lift_bias = 0.15 * exp(-(d_x^2) / (2 * (sigma/2)^2));

y_d_mod = y_d;
y_d_mod(4) = y_d(4) + lift_bias;   % Apply to swing knee
dy_d_mod = dy_d;
ddy_d_mod = ddy_d;
end


function [y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq)
% GET_ACTUALOUTPUTS Isolates the actuated degrees of freedom.
%
% Outputs:
%   y          - 4x1 vector of actual controlled outputs
%   J_y        - 4x7 Jacobian matrix of the outputs w.r.t q
%   Jdot_y_dq  - 4x1 vector representing d/dt(J_y)*dq

% Selection matrix H_0 to isolate the 4 actuated joints (last 4 states)
H_0 = [zeros(4, 3), eye(4)];

% Actual output values
y = H_0 * q;

% The output Jacobian is simply the constant selection matrix
J_y = H_0;

% Because J_y is constant, its time derivative is zero
Jdot_y_dq = zeros(4, 1);
end


function [y_d, dy_d, ddy_d] = get_DesiredOutputs(t)
% GET_DESIREDOUTPUTS Generates a smooth joint reference trajectory.
%
% Outputs:
%   y_d   - 4x1 vector of desired joint angles [rad]
%   dy_d  - 4x1 vector of desired joint velocities [rad/s]
%   ddy_d - 4x1 vector of desired joint accelerations [rad/s^2]

% Nominal standing/walking joint targets (Stance Hip, Stance Knee, Swing Hip, Swing Knee)
y0 = [0.2;  0.4; -0.2;  0.4];

% Oscillatory tracking amplitudes and tracking frequency
Amp = [0.1;  0.15;  0.1;  0.15];
omega = 2 * pi * 1.5; % 1.5 Hz walking frequency cyclic cadence

% Compute analytic trajectories, velocities, and accelerations
y_d   = y0 + Amp .* sin(omega * t);
dy_d  = Amp .* omega .* cos(omega * t);
ddy_d = -Amp .* (omega^2) .* sin(omega * t);
end