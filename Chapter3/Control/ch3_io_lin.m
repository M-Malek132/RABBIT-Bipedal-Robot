function [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, p)
%CH3_IO_LIN  Stage 4: input-output linearization of the virtual constraints.
%
%   [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, p)
%
% The outputs y have RELATIVE DEGREE 2 -- u does not appear in ydot, only in
% ydd -- so differentiating twice along xdot = f + g u gives
%
%       ydd = Lf^2 y(x) + LgLf y(x) u
%
% From ch3_outputs, ydd = Jy*ddq - (d2yd/ds2) sdot^2, and from
% ch3_control_affine, ddq = ddq_drift + ddq_in*u.  Substituting:
%
%       Lf^2 y  = Jy * ddq_drift - (d2yd/ds2) sdot^2      (ny x 1)
%       LgLf y  = Jy * ddq_in                             (ny x nu)  DECOUPLING
%
% Both are exact -- no finite differencing -- because ch3_control_affine
% returns the u-split analytically and ch3_bezier gives d2yd/ds2 analytically.
%
% If the decoupling matrix is invertible,
%
%       u = u_ff(x) + (LgLf y)^-1 mu ,     u_ff = -(LgLf y)^-1 Lf^2 y
%
% cancels the nonlinearity and leaves the clean double integrator ydd = mu.
% Every controller in stages 5-8 differs ONLY in how it picks mu.
%
% WHEN THE DECOUPLING MATRIX IS SINGULAR.  LgLfy = Jy*ddq_in loses rank where
% the choice of outputs stops being independent of the inputs.  MEASURED, what
% drives that is the PROFILE SLOPE, not the pose: dyd/ds enters Jy directly via
% Jy = H - dyd*ds_dq, and stretching alpha about alpha_0 by 100x degrades rcond
% from 5.8e-03 to 1.0e-04.
%
% Knee angle barely matters, contrary to what this comment used to claim.
% Driving the stance knee to full extension leaves rcond flat (7.0e-03 at
% q2 = 0, against 5.8e-03 at the gait's own 1.35), and a 121x121 grid over BOTH
% knees across (-pi,pi) bottoms out at rcond 2.1e-07 -- at deep double FLEXION,
% (q2,q4) = (2.13,1.87), not at lock.
%
% The pinv fallback below has never been observed to fire: over the reference
% gait rcond stays in [5.6e-03, 8.3e-03], and across 20000 random states
% spanning the full +-pi joint box the worst value is 8.0e-06, still five orders
% above the threshold.  It is genuinely defensive, not a path the pipeline
% relies on.  info.rcond reports the conditioning; callers that must not fail
% (the ODE right-hand side) should check it rather than discovering a NaN
% downstream.
%
% Inputs
%   x     : 14x1 state
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%
% Outputs
%   Lf2y  : ny x 1
%   LgLfy : ny x nu decoupling matrix
%   u_ff  : nu x 1 feedforward torque (the u that makes ydd = 0)
%   info  : struct with .y .ydot .eta .o .aux .rcond, so callers reuse this
%           single evaluation instead of recomputing the dynamics
%
% See also CH3_CONTROL_AFFINE, CH3_OUTPUTS, CH3_CTRL_PD, CH3_CTRL_CLF_QP.

[~, ~, aux]   = ch3_control_affine(x, p);
[y, ydot, o]  = ch3_outputs(x, alpha, p);

Lf2y  = o.Jy * aux.ddq_drift - o.curv;
LgLfy = o.Jy * aux.ddq_in;

rc = rcond(LgLfy);
if ~isfinite(rc) || rc < 1e-12
    % Minimum-norm fallback so an ill-conditioned pose degrades gracefully
    % instead of returning Inf/NaN into an ODE step.
    u_ff = -pinv(LgLfy) * Lf2y;
else
    u_ff = -LgLfy \ Lf2y;
end

if nargout > 3
    info = struct('y', y, 'ydot', ydot, 'eta', [y; ydot], ...
                  'o', o, 'aux', aux, 'rcond', rc);
end

end
