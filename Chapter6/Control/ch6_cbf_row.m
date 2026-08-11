function row = ch6_cbf_row(B, p)
%CH6_CBF_ROW  One barrier -> one linear inequality in u.  (6.1)/(6.3)/(6.12)
%
%   row = ch6_cbf_row(B, p)
%
% This is where the two forms Chapter 6 uses become the SAME row with different
% right-hand sides, which is the cleanest statement of how 6.1 and 6.2 relate.
%
% Both start from the relative-degree-1 function of (6.1),
%
%       h_CBF = gamma_b g + gdot,        hdot_CBF = gamma_b gdot + gddot
%
% and both demand that h_CBF decay no faster than some class-K function of
% itself. They differ only in which one:
%
%   'reciprocal'   Section 6.1: B_rec = 1/h_CBF with Bdot <= gamma/B, which
%                  after substituting B_rec is
%
%                       hdot_CBF >= -gamma h_CBF^3
%
%   'exponential'  Sections 6.2 and 6.4, (6.12):
%
%                       hdot_CBF >= -gamma h_CBF
%
% Substituting hdot_CBF = gamma_b gdot + Lf2g + LgLfg u and moving u to the
% left, BOTH are
%
%       -LgLfg u  <=  Lf2g + gamma_b gdot + kappa(h_CBF)                 (*)
%
% with kappa(h) = gamma h^3 or gamma h. So the QP row is one line and the form
% is one scalar function.
%
% ---------------------------------------------------- and it is the rb=2 ECBF
% Expanding the exponential case,
%
%       gddot >= -(gamma_b gamma) g - (gamma_b + gamma) gdot  =  -Kb eta_b
%
% with Kb = [gamma_b*gamma, gamma_b+gamma] -- which is exactly ch5_pole_gain
% for poles (gamma_b, gamma), so row.Kb is built THAT way rather than by typing
% the products out. ch6_test_barrier asserts the two agree, which is what makes
% the claim in ch6_params that "6.1.1 is the rb = 2 case of 5.2" checkable
% rather than a comment.
%
% Both poles must be real and positive: Theorem 5.2 needs Ab total negative,
% not merely Hurwitz, because Theorem 5.1's induction runs through the real
% first-order filters of (5.30). ch5_pole_gain rejects anything else.
%
% ------------------------------------------------- cubic vs linear, in practice
% Near the boundary h^3 << h, so the reciprocal form permits a MUCH faster
% approach to the constraint than the exponential one at the same gamma. That
% is the right way round for Section 6.1, which has one nominal gait and needs
% every bit of freedom it can get; it is the wrong way round for Section 6.4,
% which also has to hold a friction cone and a torque box and would rather
% start correcting early. The default in ch6_params is therefore exponential.
%
% -------------------------------------------------------- when the row is dead
% LgLfg = 0 means the row does not contain u. It is then either vacuous (rhs
% >= 0) or unsatisfiable (rhs < 0), exactly the Section 5.1 degeneracy, and the
% two are distinguished here rather than lumped together: row.live is false
% either way, row.vacuous says which. ch6_ctrl_cbf_clf_qp drops the vacuous one
% and reports the other as infeasible instead of handing quadprog a row of
% zeros with a negative right-hand side.
%
% Inputs
%   B : one element of the array from ch6_barrier / ch6_bar_lift
%   p : parameter struct (uses p.cbf)
%
% Output
%   row : struct
%           .A       1 x nu   coefficient of u, for  A u <= b
%           .b       scalar
%           .h_cbf   gamma_b g + gdot                                    (6.1)
%           .Kb      1x2      exponential-form gain, [] for reciprocal
%           .margin  h_cbf's own decrease residual, >= 0 when satisfied
%           .live    true when the row constrains u
%           .vacuous true when it does not constrain u AND is satisfied anyway
%
% See also CH6_BAR_LIFT, CH6_CTRL_CBF_CLF_QP, CH5_POLE_GAIN, CH5_ECBF_GAIN.

gb = p.cbf.gamma_b;
gm = p.cbf.gamma;

h_cbf = gb * B.g + B.gdot;                 % (6.1)

switch lower(p.cbf.form)
    case 'reciprocal'
        kappa = gm * h_cbf^3;
        Kb    = [];
    case 'exponential'
        kappa = gm * h_cbf;
        Kb    = ch5_pole_gain([gb, gm], 2);
    otherwise
        error('ch6_cbf_row:form', ...
              'Unknown p.cbf.form "%s" (expected reciprocal|exponential).', ...
              p.cbf.form);
end

rhs = B.Lf2g + gb * B.gdot + kappa;        % (*)

row = struct();
row.A       = -B.LgLfg;
row.b       = rhs;
row.h_cbf   = h_cbf;
row.Kb      = Kb;
row.live    = norm(B.LgLfg, inf) > p.cbf.min_LgLf;
row.vacuous = ~row.live && rhs >= 0;
row.label   = B.label;

% The residual the guarantee is actually about: hdot_CBF + kappa(h_CBF), which
% is >= 0 exactly when the row holds. It is filled in by the controller once u
% is known -- there is no u here -- so it starts as the u-free part.
row.margin  = NaN;

end
