function [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(x, alpha, p, unc)
%CH4_IO_LIN  I/O linearization against a CHOSEN model (nominal or true).
%
%   [...] = ch4_io_lin(x, alpha, p)        perturbed model, p.uncertainty
%   [...] = ch4_io_lin(x, alpha, p, unc)   perturbed model, given unc
%   [...] = ch4_io_lin(x, alpha, p, [])    NOMINAL model -> ch3_io_lin
%
% Identical mathematics to ch3_io_lin; the only difference is WHICH vector
% fields the Lie derivatives are taken along.  That distinction is the whole
% subject of Chapter 4, so it is made explicit in the signature rather than
% hidden in a parameter.
%
% The controller calls this with the NOMINAL model and gets the tilde
% quantities of (4.1)-(4.2):
%
%       utilde_ff = -(Lgtil Lftil y)^-1 Lftil^2 y
%       u         = utilde_ff + (Lgtil Lftil y)^-1 mu
%
% ch4_uncertainty calls it with BOTH models to form Delta1, Delta2. Nothing
% else in the chapter is allowed to ask for the true model: a controller that
% did would be cheating, and the point is to earn robustness without it.
%
% WHAT IS AND IS NOT MODEL DEPENDENT.  y, ydot, Jy and the profile curvature
% come from ch3_outputs, which is pure kinematics plus alpha -- no masses
% appear. So the outputs are IDENTICAL for both models and only
%
%       Lf2y  = Jy*ddq_drift - curv        (ddq_drift is model dependent)
%       LgLfy = Jy*ddq_in                  (ddq_in    is model dependent)
%
% differ. That is exactly why the uncertainty enters as Delta1, Delta2 acting
% on ydd in (4.3) and nowhere else.
%
% Inputs
%   x, alpha, p : as ch3_io_lin
%   unc         : perturbation; omitted -> p.uncertainty; [] -> nominal
%
% Outputs
%   Lf2y, LgLfy, u_ff, info : as ch3_io_lin, for the selected model
%
% See also CH3_IO_LIN, CH4_CONTROL_AFFINE, CH4_UNCERTAINTY.

if nargin < 4, unc = p.uncertainty; end

[~, ~, aux]  = ch4_control_affine(x, p, unc);
[y, ydot, o] = ch3_outputs(x, alpha, p);

Lf2y  = o.Jy * aux.ddq_drift - o.curv;
LgLfy = o.Jy * aux.ddq_in;

rc = rcond(LgLfy);
if ~isfinite(rc) || rc < 1e-12
    u_ff = -pinv(LgLfy) * Lf2y;
else
    u_ff = -LgLfy \ Lf2y;
end

if nargout > 3
    info = struct('y', y, 'ydot', ydot, 'eta', [y; ydot], ...
                  'o', o, 'aux', aux, 'rcond', rc);
end

end
