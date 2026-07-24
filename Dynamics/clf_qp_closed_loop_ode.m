function dxi = clf_qp_closed_loop_ode(t, xi, coeffs, theta_minus, theta_plus, p)
%CLF_QP_CLOSED_LOOP_ODE  Closed-loop swing-phase RHS under the CLF-QP law.
%
%   dxi = clf_qp_closed_loop_ode(t, xi, coeffs, theta_minus, theta_plus, p)
%
% The CLF-QP analogue of hzd_closed_loop_ode: torque comes from
% rabbit_clf_qp_controller instead of the fixed PD virtual-constraint law, and
% the returned accelerations are the QP's constrained-dynamics accelerations
% (so lambda and u are self-consistent). State xi = [x; torque_cost], and the
% last component integrates sum(tau.^2) exactly as the PD version does, so the
% same simulate/plot machinery applies.
%
% Unlike the PD closed loop this does not take a foot_ref/Baumgarte argument:
% the controller uses the acceleration-level constraint (Baumgarte off), the
% configuration the PD gait is actually tuned at (p.baumgarte_* default to 0).

    nq = p.nq;
    x  = xi(1:2*nq);
    q  = x(1:nq);
    dq = x(nq+1:end);

    [tau, info] = rabbit_clf_qp_controller(q, dq, coeffs, theta_minus, theta_plus, p);

    dxi = [dq; info.ddq; sum(tau.^2)];
end
