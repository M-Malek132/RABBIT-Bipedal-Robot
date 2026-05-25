function [u, delta_clf, delta_cbf] = rabbit_clf_cbf_controller_v2(t, x, params, target_stone_idx)
% ---------------------------------------------------------------------
% rabbit_clf_cbf_controller
% CLF-CBF controller with annulus HOCBF safety for stepping stones
% ---------------------------------------------------------------------

p  = packParameters(params);
q  = x(1:7,1);
dq = x(8:14,1);

% Robot dynamics
D = D_matrix(q,p);
C = C_vector(q,dq,p);
G = G_vector(q,p);
H_dyn = C + G;
B = input_matrix();
D_inv = inv(D);

% ---------------------------------------------------------------------
% 1. CLF Formulation
% ---------------------------------------------------------------------
[y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq);
[y_d, dy_d, ddy_d]  = get_DesiredOutputs(t);

e  = y - y_d;
de = J_y * dq - dy_d;
eta = [e; de];

% Output dynamics
Lfe_drift = Jdot_y_dq - J_y * D_inv * H_dyn - ddy_d;
Lge_ctrl  = J_y * D_inv * B;

F_drift = [de; Lfe_drift];
G_ctrl  = [zeros(length(e), size(B,2)); Lge_ctrl];

Kp_mat = 100 * eye(length(e));
Kd_mat = 20  * eye(length(e));

P = [Kp_mat + 0.5*(Kd_mat^2), 0.5*Kd_mat;
     0.5*Kd_mat,             0.5*eye(length(e))];

V   = eta' * P * eta;
LfV = 2 * eta' * P * F_drift;
LgV = 2 * eta' * P * G_ctrl;

c_clf = 10;

% ---------------------------------------------------------------------
% 2. HOCBF Formulation: Annulus Constraints
% ---------------------------------------------------------------------
stone = params.stones(target_stone_idx, :);

% Example annulus geometry
% You should tune these to match your stepping-stone geometry
R1 = 0.20;   % outer radius
R2 = 0.05;   % inner radius

% Centers of the annulus boundaries
% Example choice based on stone interval midpoint
stone_center_x = 0.5 * (stone(1) + stone(2));
O1 = [stone_center_x; 0.20];
O2 = [stone_center_x; 0.00];

% Swing foot kinematics
[~, p_sw, ~, ~, ~, ~] = rabbit_kinematics(q, p);
J_sw = J_swing(q, p);
Jdot_sw_dq = Jdotdq_swing(q, dq, p);

x_sw = p_sw(1);
z_sw = p_sw(2);

% Barrier functions
h1 = R1^2 - ((x_sw - O1(1))^2 + (z_sw - O1(2))^2);
h2 = ((x_sw - O2(1))^2 + (z_sw - O2(2))^2) - R2^2;

% Gradients wrt swing-foot position
grad_h1 = [2*(O1(1) - x_sw), 2*(O1(2) - z_sw)];
grad_h2 = [2*(x_sw - O2(1)), 2*(z_sw - O2(2))];

% First derivatives
h1_dot = grad_h1 * J_sw * dq;
h2_dot = grad_h2 * J_sw * dq;

% ---------------------------------------------------------------------
% HOCBF construction
% We use:
%   psi1 = h_dot + alpha1*h
%   psi2 = psi1_dot + alpha2*psi1 >= 0
%
% This means we need terms of the form:
%   h_ddot = Lf^2 h + LgLf h * u
%
% For implementation, we approximate the second derivative contribution
% through swing-foot acceleration dynamics.
% ---------------------------------------------------------------------

alpha1 = 10;
alpha2 = 20;

% Swing foot velocity
v_sw = J_sw * dq;   % 2x1

% Swing foot acceleration drift term:
% a_sw = Jdot*qdot + J*ddq
% ddq = D^-1 (B*u - H_dyn)
% so drift acceleration = Jdot*qdot - J*D^-1*H_dyn
a_sw_drift = Jdot_sw_dq - J_sw * D_inv * H_dyn;

% The control influence on swing-foot acceleration:
a_sw_ctrl = J_sw * D_inv * B;

% h_ddot = grad_h * a_sw + v' * Hessian(h) * v
% For squared-distance barriers, Hessian terms are constant:
%   h1 = R1^2 - ||p-O1||^2  => Hessian = -2*I
%   h2 = ||p-O2||^2 - R2^2  => Hessian = +2*I
%
% Therefore:
h1_ddot_drift = grad_h1 * a_sw_drift - 2 * (v_sw' * v_sw);
h2_ddot_drift = grad_h2 * a_sw_drift + 2 * (v_sw' * v_sw);

h1_ddot_ctrl = grad_h1 * a_sw_ctrl;
h2_ddot_ctrl = grad_h2 * a_sw_ctrl;

% Since h_ddot = h_ddot_drift + h_ddot_ctrl*u
% and psi2 >= 0 gives:
%   h_ddot + alpha2*h_dot + alpha1*h >= 0
%
% So constraint becomes:
%   (h_ddot_ctrl) u >= -h_ddot_drift - alpha2*h_dot - alpha1*h

A_cbf = -[h1_ddot_ctrl;
          h2_ddot_ctrl];

b_cbf = [h1_ddot_drift + alpha2*h1_dot + alpha1*h1;
         h2_ddot_drift + alpha2*h2_dot + alpha1*h2];

% ---------------------------------------------------------------------
% 3. Quadratic Program
% ---------------------------------------------------------------------
% Decision variables:
% z = [u1; u2; u3; u4; delta_clf; delta_cbf]
%
% Objective:
% minimize ||u||^2 + p_clf*delta_clf^2 + p_cbf*delta_cbf^2

weight_u = 1;
weight_delta_clf = 1000;
weight_delta_cbf = 1000;

H_qp = diag([weight_u*ones(4,1); weight_delta_clf; weight_delta_cbf]);
f_qp = zeros(6,1);

% CLF constraint:
% LfV + LgV*u <= -c_clf*V + delta_clf
% => [LgV, -1, 0] * z <= -LfV - c_clf*V
A_clf = [LgV, -1, 0];
b_clf = -LfV - c_clf * V;

% CBF constraints:
% A_cbf*u <= b_cbf
% => [A_cbf, 0, 0] * z <= b_cbf
A_cbf_qp = [A_cbf, zeros(size(A_cbf,1), 2)];

A_ineq = [A_clf;
          A_cbf_qp];

b_ineq = [b_clf;
          b_cbf];

% Input bounds
u_max = 150;
lb = [-u_max*ones(4,1); 0; 0];
ub = [ u_max*ones(4,1); Inf; Inf];

% Solve QP
options = optimoptions('quadprog', 'Display', 'off');
[z_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);

if exitflag == 1
    u = z_opt(1:4);
    delta_clf = z_opt(5);
    delta_cbf = z_opt(6);
else
    % -------------------------------------------------------------
    % Fallback: PD control if QP fails
    % -------------------------------------------------------------
    warning('QP infeasible at t=%.3f. Falling back to PD control.', t);

    Kp_pd = 150;
    Kd_pd = 20;

    error = y_d - y;
    d_error = dy_d - (J_y * dq);

    u = Kp_pd * error + Kd_pd * d_error;
    u = max(min(u, u_max), -u_max);

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
