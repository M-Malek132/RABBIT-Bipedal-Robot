function [u, ci] = ch6_control(t, x, alpha, p)
%CH6_CONTROL  Single dispatch point for the Chapter-6 feedback law.
%
%   [u, ci] = ch6_control(t, x, alpha, p)
%
% Same role as ch3_control: every consumer of a torque routes through here, so
% p.controller changes the whole pipeline with no other edit. Two differences
% from Chapter 3, both forced by the chapter:
%
%   1. IT TAKES TIME. Chapter 3's feedback is time-invariant -- that is the
%      point of a phase variable. Section 6.2's stepping stones MOVE, so the
%      barrier depends on t explicitly and the signature has to carry it. t is
%      time since the START OF THE CURRENT STEP, not global time: the stone
%      freezes on impact and the next stone starts its own clock.
%
%   2. IT TAKES alpha PER CALL, and Section 6.4 changes alpha per STEP. The
%      gait library hands a freshly interpolated alpha to each step; within a
%      step it is constant, which is what keeps the virtual constraints
%      well-defined and the outputs continuous.
%
% Controllers:
%   'ff'          u = u_ff.  No feedback, no barrier.
%   'iolin_pd'    Chapter 3 stage 5.
%   'clfqp'       Chapter 3 stage 7.
%   'clfqp_con'   Chapter 3 stage 8 -- the no-CBF baseline, controller I of
%                 (6.27) when the gait library is on.
%   'cbf_clf_qp'  Chapter 6, (6.12)/(6.26).
%
% The Chapter-3 names are delegated to ch3_control rather than reimplemented,
% so a change there cannot make the Table 6.1 baselines drift away from the
% controllers they are baselines for.
%
% Inputs
%   t     : time since the start of the current step [s]
%   x     : 14x1 state
%   alpha : ny x n_ctrl Bezier coefficients for THIS step
%   p     : parameter struct
%
% Outputs
%   u  : nu x 1 joint torque
%   ci : control info -- everything ch3_control returns, plus
%          .B      the barrier struct array
%          .qp     the Chapter-6 QP report (see ch6_ctrl_cbf_clf_qp)
%          .geom   problem geometry for plots
%          .h      1 x nb barrier values
%
% See also CH6_CTRL_CBF_CLF_QP, CH6_BARRIER, CH3_CONTROL.

if ~strcmpi(p.controller, 'cbf_clf_qp')
    % Chapter-3 law, no barrier rows. Still reports the barrier VALUES so a
    % baseline run can be plotted on the same axes as a Chapter-6 run and the
    % violation is visible rather than inferred.
    [u, ci] = ch3_control(x, alpha, p);

    if nargout > 1
        [B, geom] = ch6_barrier(x, ci.aux, p, t);
        ci.B    = B;
        ci.geom = geom;
        ci.h    = arrayfun(@(b) b.g, B);
        ci.qp   = struct('feasible', ci.qp_feasible, 'delta', ci.delta, ...
                         'h', ci.h, 'margin', nan(1, numel(B)), ...
                         'cbf_active', false(1, numel(B)), 'n_dead', 0, ...
                         'viol', 0);
    end
    return;
end

[Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, p);
[B, geom] = ch6_barrier(x, info.aux, p, t);

% WHAT THE QP'S COST PULLS TOWARDS. See ch6_ctrl_cbf_clf_qp for the measurement
% behind the default: the literal min-norm reading of (6.26) destabilises this
% gait in four steps WITH NO BARRIERS AT ALL, so it is not the barrier that
% needs the margin, it is the cost.
switch lower(p.cbf.reference)
    case 'ff'
        u_ref = u_ff;                                  % mu_ref = 0, literal (6.26)
    case 'pd'
        [~, u_ref] = ch3_ctrl_pd(Lf2y, LgLfy, u_ff, info, p);
    otherwise
        error('ch6_control:reference', ...
              'Unknown p.cbf.reference "%s" (expected ff|pd).', p.cbf.reference);
end

[u, qp] = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, B, p, u_ref);

if nargout > 1
    ci = struct('y', info.y, 'ydot', info.ydot, 'eta', info.eta, ...
                'mu', qp.mu, 'u_ff', u_ff, 'Lf2y', Lf2y, 'LgLfy', LgLfy, ...
                'aux', info.aux, 'o', info.o, 'rcond', info.rcond, ...
                'V', qp.V, 'delta', qp.delta, 'qp_feasible', qp.feasible, ...
                'B', B, 'qp', qp, 'geom', geom, 'h', qp.h);
end

end
