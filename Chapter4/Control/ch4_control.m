function [u, xidot, ci] = ch4_control(x, xi, alpha, p)
%CH4_CONTROL  Single dispatch point for every Chapter-4 feedback law.
%
%   [u, xidot, ci] = ch4_control(x, xi, alpha, p)
%
% Plays the same role as ch3_control -- one place that every torque consumer
% routes through -- with two differences that define the chapter.
%
% FIRST: THE CONTROLLER IS BUILT FROM THE NOMINAL MODEL, ALWAYS.  The single
% call to ch4_io_lin below passes unc = [], so every law here sees ftil, gtil.
% The true model enters only in ch4_ode_rhs, which integrates the plant. That
% separation is enforced structurally rather than by convention: there is no
% path from this function to the perturbed dynamics.
%
% SECOND: CONTROLLERS MAY HAVE STATE.  The L1 laws carry a predictor, two
% parameter estimates and a filter (ch4_l1_state), so the signature takes xi
% and returns xidot. Stateless laws return xidot = [] and ignore xi, so the
% caller does not need to know which is which.
%
% THE MENU
%   'ff' 'iolin_pd' 'clfqp' 'clfqp_con'   Chapter 3, unchanged. Under nonzero
%                                         uncertainty these are the BASELINES
%                                         (controller A in Sections 4.1.4 and
%                                         4.2.4) -- they do not know the model
%                                         is wrong and cannot react to it.
%   'rclfqp'                              robust CLF-QP, eq (4.12)
%   'rclfqp_con'                          robust CLF-QP + constraints, (4.13)
%   'l1'                                  L1 adaptive + CLF-QP, Section 4.2.2
%   'l1_con'                              the same + torque saturation, 4.2.3
%
% Inputs
%   x     : 14x1 state
%   xi    : controller state (5ny x 1 for L1, [] otherwise)
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%
% Outputs
%   u     : nu x 1 joint torque
%   xidot : controller state derivative, [] for stateless laws
%   ci    : control info, uniform across all laws --
%             .y .ydot .eta .mu .u_ff .Lf2y .LgLfy .aux .o .rcond
%             .V .delta .qp_feasible
%             .qp   robust-QP detail, or []
%             .l1   L1 detail (theta_hat, eta_tilde, mu1, mu2, ...), or []
%
%           NOTE .aux is the NOMINAL model's dynamics data. Do not use it to
%           report ground reaction forces -- those must come from the true
%           model. ch4_forces does that correctly; see its header.
%
% See also CH3_CONTROL, CH4_IO_LIN, CH4_CTRL_RCLF_QP, CH4_CTRL_L1.

[Lf2y, LgLfy, u_ff, info] = ch4_io_lin(x, alpha, p, []);   % NOMINAL, always

V = NaN; delta = 0; qp_ok = true; qp = []; l1 = [];
xidot = [];

switch lower(p.controller)

    case 'ff'
        mu = zeros(p.ny, 1);
        u  = u_ff;

    case 'iolin_pd'
        [mu, u] = ch3_ctrl_pd(Lf2y, LgLfy, u_ff, info, p);

    case 'clfqp'
        [mu, u, qp] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p, false);
        V = qp.V;

    case 'clfqp_con'
        [mu, u, qp] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p, true);
        V = qp.V; delta = qp.delta; qp_ok = qp.feasible;

    case 'rclfqp'
        [mu, u, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, p, false);
        V = qp.V; qp_ok = qp.feasible;

    case 'rclfqp_con'
        [mu, u, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, p, true);
        V = qp.V; delta = qp.delta; qp_ok = qp.feasible;

    case 'l1'
        [mu, u, xidot, l1] = ch4_ctrl_l1(Lf2y, LgLfy, u_ff, info, xi, p, false);
        V = l1.V; qp_ok = l1.qp_feasible;

    case 'l1_con'
        [mu, u, xidot, l1] = ch4_ctrl_l1(Lf2y, LgLfy, u_ff, info, xi, p, true);
        V = l1.V; delta = l1.delta; qp_ok = l1.qp_feasible;

    otherwise
        error('ch4_control:controller', ...
              ['Unknown p.controller "%s" (expected ff|iolin_pd|clfqp|' ...
               'clfqp_con|rclfqp|rclfqp_con|l1|l1_con).'], p.controller);
end

if nargout > 2
    ci = struct('y', info.y, 'ydot', info.ydot, 'eta', info.eta, ...
                'mu', mu, 'u_ff', u_ff, 'Lf2y', Lf2y, 'LgLfy', LgLfy, ...
                'aux', info.aux, 'o', info.o, 'rcond', info.rcond, ...
                'V', V, 'delta', delta, 'qp_feasible', qp_ok, ...
                'qp', qp, 'l1', l1);
end

end
