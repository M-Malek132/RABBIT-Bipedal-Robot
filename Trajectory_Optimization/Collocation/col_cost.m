function J = col_cost(z, p)
%COL_COST  Torque-squared per unit step length (thesis Eq. 3.5 / Westervelt 6.43).
%
%       J = (1/L_step) * integral_0^T ||u(t)||^2 dt
%
% The integral uses the SAME Hermite-Simpson quadrature as the dynamics
% transcription (midpoint torque = node average), so cost and defects are
% mutually consistent. Step length is clamped and a short-step quadratic
% penalty is added, mirroring the shooting cost (hzd_cost) so the two solvers
% optimize the same objective.

    [X, U, ~, ~, T] = col_unpack(z, p);
    N = p.N;
    h = T/(N-1);

    g = sum(U.^2, 1);                 % ||u_k||^2 at each node, 1 x N
    Ig = 0;
    for k = 1:N-1
        umid = 0.5*(U(:,k) + U(:,k+1));
        gmid = sum(umid.^2);
        Ig = Ig + (h/6)*(g(k) + 4*gmid + g(k+1));
    end

    qN     = X(1:p.nq, N);
    sw     = P_sw(qN);   st = P_st(qN);
    L_step = sw(1) - st(1);

    L_eff = max(L_step, p.min_step_length);
    J = Ig / L_eff;
    J = J + p.short_step_weight * max(0, p.min_step_length - L_step)^2;
end
