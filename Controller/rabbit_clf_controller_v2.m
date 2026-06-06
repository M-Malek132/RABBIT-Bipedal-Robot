% -------------------------------------------------------------------------
% Robust CLF-QP Controller
% -------------------------------------------------------------------------
function [u, delta_clf] = rabbit_clf_controller_v2(t, x, params)
% Persistent variables to track trajectory anchoring across hybrid steps
persistent x_initial last_step_idx

global CURRENT_STEP;

% Reset the anchor point if we have transitioned to a new step
if isempty(x_initial) || isempty(last_step_idx) || (CURRENT_STEP ~= last_step_idx)
    x_initial = x;
    last_step_idx = CURRENT_STEP;
end

%% 1. Extract Dynamics
p  = packParameters(params);
q  = x(1:7,1);
dq = x(8:14,1);

D = D_matrix(q,p);
C = C_vector(q,dq,p);
G = G_vector(q,p);
H_total = C + G;
B = input_matrix();

%% 2. Tracking Error (Anchored to x_initial)
[y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq);
[y_d, dy_d, ddy_d]  = get_DesiredOutputs(t, x_initial, p);

e   = y - y_d;
de  = J_y * dq - dy_d;
eta = [e; de];

%% 3. Lie Derivatives
D_inv_H = D \ H_total;
D_inv_B = D \ B;

% Drift and Control influence on the output acceleration
Lf_drift = Jdot_y_dq - J_y * D_inv_H - ddy_d;
Lg_ctrl  = J_y * D_inv_B;

%% 4. Control Lyapunov Function (CLF)
n_y = length(e);
Kp_val = 150; % Increased for tighter tracking
Kd_val = 25;
Kp_mat = Kp_val * eye(n_y);
Kd_mat = Kd_val * eye(n_y);

% Lyapunov matrix P (Solution to F'P + PF = -Q)
P = [Kp_mat + 0.5*(Kd_mat^2), 0.5*Kd_mat;
     0.5*Kd_mat,             0.5*eye(n_y)];

% State space error dynamics: d_eta = F_drift + G_ctrl * u
F_err = [de; Lf_drift];
G_err = [zeros(n_y, size(B,2)); Lg_ctrl];

V   = eta' * P * eta;
LfV = 2 * eta' * P * F_err;
LgV = 2 * eta' * P * G_err;

% Convergence rate (soft start for first 0.1s of the step)
c_clf = 5.0; 
if t < 0.1, c_clf = 1.0; end

%% 5. QP Formulation
n_u = size(B,2);
n_z = n_u + 1; % [u; delta_clf]

weight_u = 0.01;
weight_delta_clf = 2000; % High penalty on CLF relaxation

H_qp = diag([weight_u * ones(1, n_u), weight_delta_clf]);
H_qp = H_qp + eye(n_z) * 1e-9; % Numerical regularization
f_qp = zeros(n_z, 1);

% CLF Inequality: LfV + LgV*u <= -c_clf * V + delta_clf
A_ineq = [LgV, -1];
b_ineq = -LfV - c_clf * V;

% Actuator Constraints
u_max = 150; 
lb = [-u_max * ones(n_u, 1); 0];
ub = [ u_max * ones(n_u, 1); Inf];

options = optimoptions('quadprog', ...
    'Algorithm', 'interior-point-convex', ...
    'Display', 'off', ...
    'ConstraintTolerance', 1e-4, ...
    'OptimalityTolerance', 1e-4);

[z_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);

%% 6. Output and Fallback
if exitflag == 1 || exitflag == 0
    u = z_opt(1:n_u);
    delta_clf = z_opt(n_z);
else
    % Safe Fallback: Saturated PD Control
    u_pd = -(Kp_val*0.1 * e + Kd_val*0.2 * de);
    u = max(min(u_pd, u_max), -u_max);
    delta_clf = 0;
end
end
