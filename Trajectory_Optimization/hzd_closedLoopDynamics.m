function [dx, u] = hzd_closedLoopDynamics(t, x, model, opt, simOpt)
%HZD_CLOSEDLOOPDYNAMICS  Swing-phase ODE rhs with HZD controller.
%
%  Reads simOpt.CP  (control-point matrix, (n+1)×4)
%  Calls: D_matrix, C_vector, G_vector, input_matrix  from your Dynamics/

nq     = model.nq;
params = model.params;

q  = x(1:nq);
dq = x(nq+1:end);

CP = simOpt.CP;    % (n+1)×4 — replaces simOpt.alpha
Kp = simOpt.Kp;
Kd = simOpt.Kd;

% ---- Dynamics matrices ----------------------------------------
D = D_matrix(q, params);
C = C_vector(q, dq, params);
G = G_vector(q, params);
B = input_matrix(q, params);   % 7×4

% ---- Virtual constraints -------------------------------------
[y, dy, Jy] = hzd_virtualConstraints(q, dq, CP, model, opt);

% ---- Numerical  dot(Jy)*dq  ----------------------------------
eps_fd = 1e-7;
q_p    = q + eps_fd * dq;
[~, ~, Jy_p] = hzd_virtualConstraints(q_p, dq, CP, model, opt);
Jydot_dq = ((Jy_p - Jy) / eps_fd) * dq;

% ---- Input-output linearisation  ddy = v ---------------------
v     = -Kd * dy - Kp * y;
A_mat = Jy  * (D \ B);                          % 4×4
b_vec = v - Jydot_dq + Jy * (D \ (C*dq + G));  % 4×1

u = A_mat \ b_vec;
u = max(min(u, opt.uMax), opt.uMin);

% ---- State derivative ----------------------------------------
ddq = D \ (B*u - C*dq - G);
dx  = [dq; ddq];
end
