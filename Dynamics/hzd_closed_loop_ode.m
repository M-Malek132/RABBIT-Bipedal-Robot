function dxi = hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p)
    nq = p.nq;
    x  = xi(1:2*nq);
    q  = x(1:nq);
    dq = x(nq+1:end);

    theta     = theta_of_q(q);
    s         = max(0, min(1, (theta - theta_minus) / (theta_plus - theta_minus)));
    ds_dtheta = 1/(theta_plus - theta_minus);
    dtheta_dt = dtheta_dq_of(q) * dq;
    dtheta_dt = max(-100, min(100, dtheta_dt));   % debug safety net

    act_idx = 4:7;
    y  = q(act_idx);
    dy = dq(act_idx);

    yd  = zeros(p.nu,1);
    dyd = zeros(p.nu,1);
    for i = 1:p.nu
        [b, db] = bspline_eval(coeffs(i,:), s, p);
        yd(i)  = b;
        dyd(i) = db * ds_dtheta * dtheta_dt;
    end

    e  = max(-10,  min(10,  y  - yd));   % debug safety net
    de = max(-100, min(100, dy - dyd));  % debug safety net

    tau = -p.Kp.*e - p.Kd.*de;
    tau = max(-1000, min(1000, tau));    % debug safety net

    ddq = rabbit_constrained_dynamics(q, dq, tau);
    ddq = max(-1000, min(1000, ddq));    % debug safety net

    dxi = [dq; ddq; sum(tau.^2)];
end