function Jh0 = hzd_jacobian_h0(q, model)
%HZD_JACOBIAN_H0  Jacobian of h0(q) = q(4:7) w.r.t. q  (ny x nq).
%
%  h0 = q(4:7)  =>  dh0/dq = [0_{4x3}, I_{4x4}]

nq = model.nq;   % 7
ny = model.nu;   % 4

Jh0 = zeros(ny, nq);
Jh0(:, 4:7) = eye(ny);
end
