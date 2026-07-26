function [u, ci] = ch3_control(x, alpha, p)
%CH3_CONTROL  Single dispatch point for the swing-phase feedback law.
%
%   [u, ci] = ch3_control(x, alpha, p)
%
% Every consumer of a torque -- the ODE right-hand side, the collocation
% dynamics, the cost, the force recovery in the report -- calls THIS function,
% so p.controller changes the entire pipeline consistently with no other edit.
% That single-point-of-truth property is what makes "solve under PD, then
% re-check under the constrained CLF-QP" a one-line experiment.
%
%   'ff'        u = u_ff only.  No output feedback: the exact torque that
%               holds ydd = 0. This is what the collocation solve uses, since
%               there the virtual constraints are enforced as equality
%               constraints at every node and eta is identically zero -- any
%               feedback term would contribute nothing but numerical noise.
%   'iolin_pd'  stage 5.
%   'clfqp'     stage 7, unconstrained.
%   'clfqp_con' stage 8, with the torque box and whichever of friction / GRF
%               are enabled in p.limits.enable.
%
% Inputs
%   x     : 14x1 state
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%
% Outputs
%   u  : nu x 1 joint torque
%   ci : control info -- .y .ydot .eta .mu .u_ff .Lf2y .LgLfy .aux .rcond
%        plus .V .delta .qp_feasible for the CLF paths. Callers use .aux to
%        recover the ground reaction force without redoing the KKT solve.
%
% See also CH3_IO_LIN, CH3_CTRL_PD, CH3_CTRL_CLF_QP.

[Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, p);

V = NaN; delta = 0; qp_ok = true;

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
        V     = qp.V;
        delta = qp.delta;
        qp_ok = qp.feasible;

    otherwise
        error('ch3_control:controller', ...
              ['Unknown p.controller "%s" ' ...
               '(expected ff|iolin_pd|clfqp|clfqp_con).'], p.controller);
end

if nargout > 1
    ci = struct('y', info.y, 'ydot', info.ydot, 'eta', info.eta, ...
                'mu', mu, 'u_ff', u_ff, 'Lf2y', Lf2y, 'LgLfy', LgLfy, ...
                'aux', info.aux, 'o', info.o, 'rcond', info.rcond, ...
                'V', V, 'delta', delta, 'qp_feasible', qp_ok);
end

end
