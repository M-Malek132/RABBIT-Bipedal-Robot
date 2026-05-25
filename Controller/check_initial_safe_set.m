function info = check_initial_safe_set(t0, x0, params, target_stone_idx)
p  = packParameters(params);
q  = x0(1:7,1);
dq = x0(8:14,1);

% Dynamics
D = D_matrix(q,p);
C = C_vector(q,dq,p);
G = G_vector(q,p);
H_dyn = C + G;
B = input_matrix();
D_inv = inv(D);

% Stone geometry (same as controller)
stone = params.stones(target_stone_idx, :);
R1 = 0.20;
R2 = 0.05;

stone_center_x = 0.5 * (stone(1) + stone(2));
O1 = [stone_center_x; 0.20];
O2 = [stone_center_x; 0.00];

% Swing-foot kinematics
[~, p_sw, ~, ~, ~, ~] = rabbit_kinematics(q, p);
J_sw = J_swing(q, p);
Jdot_sw_dq = Jdotdq_swing(q, dq, p);

x_sw = p_sw(1);
z_sw = p_sw(2);

% Barriers
h1 = R1^2 - ((x_sw - O1(1))^2 + (z_sw - O1(2))^2);
h2 = ((x_sw - O2(1))^2 + (z_sw - O2(2))^2) - R2^2;

grad_h1 = [2*(O1(1) - x_sw), 2*(O1(2) - z_sw)];
grad_h2 = [2*(x_sw - O2(1)), 2*(z_sw - O2(2))];

v_sw = J_sw * dq;   % 2x1

h1_dot = grad_h1 * v_sw;
h2_dot = grad_h2 * v_sw;

alpha1 = 10;
alpha2 = 20;

psi1_1 = h1_dot + alpha1*h1;
psi1_2 = h2_dot + alpha1*h2;

% accel decomposition
a_sw_drift = Jdot_sw_dq - J_sw * D_inv * H_dyn;
a_sw_ctrl  = J_sw * D_inv * B;

h1_ddot_drift = grad_h1 * a_sw_drift - 2 * (v_sw' * v_sw);
h2_ddot_drift = grad_h2 * a_sw_drift + 2 * (v_sw' * v_sw);

h1_ddot_ctrl = grad_h1 * a_sw_ctrl; % 1x4
h2_ddot_ctrl = grad_h2 * a_sw_ctrl; % 1x4

% Your QP uses:  A_cbf * u <= b_cbf  with A_cbf = -hddot_ctrl
A_cbf = -[h1_ddot_ctrl;
h2_ddot_ctrl];     % 2x4

b_cbf = [h1_ddot_drift + alpha2*h1_dot + alpha1*h1;
h2_ddot_drift + alpha2*h2_dot + alpha1*h2];  % 2x1

% Check feasibility with input bounds:
u_max = inf; % change this to your real saturation if you have it
lb_u = -u_max*ones(4,1);
ub_u =  u_max*ones(4,1);

% Solve a feasibility LP/QP: find u such that A_cbf u <= b_cbf
% Use linprog if you have finite bounds; with inf bounds, this is unbounded but still feasible-checkable.
H = eye(4); f = zeros(4,1);
opts = optimoptions('quadprog','Display','off');
[~,~,exitflag] = quadprog(H,f,A_cbf,b_cbf,[],[],lb_u,ub_u,[],opts);

info = struct();
info.t0 = t0;
info.x_sw = x_sw; info.z_sw = z_sw;
info.h1 = h1; info.h2 = h2;
info.h1_dot = h1_dot; info.h2_dot = h2_dot;
info.psi1_1 = psi1_1; info.psi1_2 = psi1_2;
info.A_cbf = A_cbf; info.b_cbf = b_cbf;
info.cbf_feasible_with_bounds = (exitflag == 1);

% disp(info)
end