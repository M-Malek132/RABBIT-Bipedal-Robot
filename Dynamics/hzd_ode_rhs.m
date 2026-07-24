function dxi = hzd_ode_rhs(t, xi, coeffs, theta_minus, theta_plus, p, foot_ref)
%HZD_ODE_RHS  Swing-phase closed-loop RHS, dispatched on the chosen controller.
%
%   dxi = hzd_ode_rhs(t, xi, coeffs, theta_minus, theta_plus, p [, foot_ref])
%
% Single source of truth for WHICH feedback law the single-support simulation
% (and therefore the gait optimizer) runs. Selected by p.controller:
%
%   'pd'    (default)  hzd_closed_loop_ode      -- fixed-gain virtual-constraint
%                                                  PD law (uses foot_ref for the
%                                                  optional Baumgarte term).
%   'clfqp'            clf_qp_closed_loop_ode   -- input-output linearization +
%                                                  control-Lyapunov-function QP
%                                                  (Chapter 3). Respects the
%                                                  torque box / friction cone
%                                                  that the PD law cannot.
%
% Any unset/unrecognized p.controller falls back to 'pd', so existing callers
% and saved p-structs that predate this field are unaffected.

    if nargin < 7, foot_ref = []; end

    if isfield(p, 'controller') && strcmpi(p.controller, 'clfqp')
        dxi = clf_qp_closed_loop_ode(t, xi, coeffs, theta_minus, theta_plus, p);
    else
        dxi = hzd_closed_loop_ode(t, xi, coeffs, theta_minus, theta_plus, p, foot_ref);
    end
end
