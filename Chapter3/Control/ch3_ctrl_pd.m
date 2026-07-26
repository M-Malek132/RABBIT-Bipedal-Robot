function [mu, u] = ch3_ctrl_pd(Lf2y, LgLfy, u_ff, info, p)
%CH3_CTRL_PD  Stage 5: the baseline -- PD on the linearized system.
%
%   [mu, u] = ch3_ctrl_pd(Lf2y, LgLfy, u_ff, info, p)
%
% Stage 4 left a clean double integrator ydd = mu.  The simplest thing to put
% there is a PD law, scaled by eps:
%
%       mu = -(1/eps^2) Kp y - (1/eps) Kd ydot
%
% The closed loop is etadot = A_eps eta, and it is stable iff
%
%       [   0     I  ]
%       [ -Kp   -Kd  ]      is Hurwitz.
%
% eps does NOT change stability -- it changes the RATE.  The eps-scaling makes
% the closed-loop eigenvalues scale like 1/eps, and that matters because the
% impact map EXPANDS eta at every step: convergence has to be fast enough,
% within one step, to beat that expansion.  Too large an eps and the gait is
% locally unstable at the step-to-step level even though the continuous phase
% looks perfectly well behaved.
%
% What this law CANNOT do -- and the reason stages 6-8 exist -- is respect a
% torque limit.  u is whatever comes out of the inversion; if that exceeds the
% actuators, the real robot silently saturates and the guarantee is void.
%
% Inputs
%   Lf2y, LgLfy, u_ff, info : from ch3_io_lin
%   p                       : parameter struct (uses p.eps, p.Kp, p.Kd)
%
% Outputs
%   mu : ny x 1 virtual input
%   u  : nu x 1 joint torque, u = u_ff + (LgLfy)^-1 mu
%
% See also CH3_IO_LIN, CH3_CTRL_CLF_QP.

y    = info.y;
ydot = info.ydot;

mu = -(1/p.eps^2) * (p.Kp * y) - (1/p.eps) * (p.Kd * ydot);

if nargout > 1
    u = u_ff + solve_decoupling(LgLfy, mu, info);
end

end

function v = solve_decoupling(A, mu, info)
if isfield(info,'rcond') && (~isfinite(info.rcond) || info.rcond < 1e-12)
    v = pinv(A) * mu;
else
    v = A \ mu;
end
end
