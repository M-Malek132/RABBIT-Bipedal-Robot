function dxi = hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p,foot_ref)
    % Closed-loop swing-phase RHS. The control law and constrained dynamics
    % now live in hzd_control_and_dynamics (so tau/lambda can be recomputed
    % post-hoc); this wrapper only appends the running torque-cost integrand.
    %
    % foot_ref (optional 2x1): world-frame stance-foot pin location; enables
    % Baumgarte stabilization in the constrained dynamics. Omit for the old
    % (acceleration-level-only) behaviour.
    if nargin < 7, foot_ref = []; end
    nq = p.nq;
    x  = xi(1:2*nq);
    q  = x(1:nq);
    dq = x(nq+1:end);

    [ddq, tau] = hzd_control_and_dynamics(q, dq, coeffs, theta_minus, theta_plus, p, foot_ref);

    dxi = [dq; ddq; sum(tau.^2)];
end
