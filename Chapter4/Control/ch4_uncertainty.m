function D = ch4_uncertainty(x, alpha, p, unc)
%CH4_UNCERTAINTY  The induced model uncertainty Delta1, Delta2 of eq (4.4).
%
%   D = ch4_uncertainty(x, alpha, p)
%   D = ch4_uncertainty(x, alpha, p, unc)
%
% Applying the NOMINAL pre-control (4.1) to the TRUE plant does not give the
% double integrator ydd = mu that Chapter 3 relied on. It gives (4.3):
%
%       ydd = mu + Delta1 + Delta2 mu
%
% with, from (4.4),
%
%       Delta1 = Lf^2 y - LgLf y (Lgtil Lftil y)^-1 Lftil^2 y
%       Delta2 = LgLf y (Lgtil Lftil y)^-1 - I
%
% where untilded quantities are the TRUE model and tilded ones the NOMINAL.
% This function is DIAGNOSTIC: it needs the true model, which no controller in
% this chapter is allowed to see. Its job is to measure what the controllers
% must survive.
%
% A CHEAPER AND MORE MEANINGFUL FORM FOR Delta1.  Since
% utilde_ff = -(Lgtil Lftil y)^-1 Lftil^2 y, the definition collapses to
%
%       Delta1 = Lf^2 y + LgLf y * utilde_ff
%
% which is literally "the ydd the true robot produces when you feed it the
% nominal feedforward". Chapter 3 designed utilde_ff precisely to make that
% zero; Delta1 is how badly that fails. This is the form computed below -- one
% matrix solve instead of two, and it makes the interpretation obvious.
%
% Remark 4.1 says Delta1 = Delta2 = 0 when the models agree. Here that is not
% an approximation to within a tolerance: with unc trivial, utilde_ff is the
% exact solution of LgLfy*u = -Lf2y, so both expressions vanish to round-off.
% ch4_test_model asserts it.
%
% WHY BOTH TERMS MATTER, DIFFERENTLY.  Delta1 is an additive disturbance: it
% destroys the equilibrium, so the closed loop (4.5) has no fixed point and
% tracking errors no longer go to zero. Delta2 multiplies the control: it can
% rotate or shrink the control's effect, and if it is large enough to reverse a
% direction, the controller actively destabilizes what it is trying to fix.
% That is why the robust QP's feasibility condition is a bound on Delta2 alone
% -- see ch4_ctrl_rclf_qp.
%
% Inputs
%   x, alpha, p : state, virtual constraint coefficients, parameters
%   unc         : perturbation defining the TRUE model; default p.uncertainty
%
% Outputs
%   D : struct with
%         .Delta1     ny x 1
%         .Delta2     ny x ny
%         .n1         ||Delta1||_2
%         .n2         ||Delta2||_2   (induced 2-norm, matches p.rclf 'matrix')
%         .n2_scalar  the scalar-model magnitude, max_i |sum_j Delta2(i,j)| /
%                     ... see below; matches p.rclf 'scalar'
%         .u_ff_nom   the nominal feedforward that produced Delta1
%         .Lf2y_true .LgLfy_true .Lf2y_nom .LgLfy_nom
%
% See also CH4_DELTA_BOUNDS, CH4_CTRL_RCLF_QP, CH4_IO_LIN.

if nargin < 4, unc = p.uncertainty; end

[Lf2y_n, LgLfy_n, u_ff_n] = ch4_io_lin(x, alpha, p, []);      % NOMINAL
[Lf2y_t, LgLfy_t]         = ch4_io_lin(x, alpha, p, unc);     % TRUE

Delta1 = Lf2y_t + LgLfy_t * u_ff_n;

% Delta2 = LgLfy_true * inv(LgLfy_nom) - I, solved as a right division so the
% nominal decoupling matrix is never explicitly inverted.
Delta2 = (LgLfy_t / LgLfy_n) - eye(p.ny);

D = struct();
D.Delta1     = Delta1;
D.Delta2     = Delta2;
D.n1         = norm(Delta1, 2);
D.n2         = norm(Delta2, 2);

% The 'scalar' uncertainty model of ch4_ctrl_rclf_qp assumes Delta2 = d2*I.
% The honest scalar summary of a general Delta2 under that model is the gain it
% applies along the direction the CLF actually cares about, which is not known
% here; the isotropic part is the model-free stand-in, and the gap between this
% and n2 is how anisotropic the real uncertainty is.
D.n2_scalar  = abs(trace(Delta2)) / p.ny;

D.u_ff_nom   = u_ff_n;
D.Lf2y_true  = Lf2y_t;
D.LgLfy_true = LgLfy_t;
D.Lf2y_nom   = Lf2y_n;
D.LgLfy_nom  = LgLfy_n;

end
