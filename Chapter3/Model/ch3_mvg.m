function [Mq, Vv, Gv] = ch3_mvg(q, dq, p)
%CH3_MVG  M(q), V(q,dq) and G(q) of the model p describes.
%
%   [Mq, Vv, Gv] = ch3_mvg(q, dq, p)
%   Mq           = ch3_mvg(q, [], p)         (Vv = [] when dq is empty)
%
% Every Chapter-3 use of the equations of motion goes through here: the
% continuous phase (ch3_control_affine), the reset map (ch3_impact) and the
% zero dynamics (ch3_zd_point). With p.model_blend empty -- the default --
% this is M.m / V.m / G.m and nothing else, so no existing result changes.
%
% p.model_blend = lambda in [0, 1] selects the robot whose link masses and
% inertias are (1 - lambda) times those of the 30 kg model this repository
% used before 2026-09-02 plus lambda times today's 74 kg ones. The geometry of
% the two is identical, and M, V and G are LINEAR in the masses and inertias
% (M is a sum of traces linear in each link's pseudo-inertia, V is built from
% M, G is linear in the masses), so blending the two generated models IS the
% model with blended parameters -- exactly, not as an approximation.
%
%   lambda = 1   today's robot (takes the plain M.m path)
%   lambda = 0   the old one: M_m30 / V_m30 / G_m30, restored from git
%
% ch3_model_homotopy walks lambda from 0 to 1 with the torque box fixed, to
% find where a gait inside the declared limit stops existing.
%
% Inputs
%   q  : 7x1 configuration
%   dq : 7x1 velocity, or [] when only M (and G) are needed
%   p  : parameter struct (uses p.model_blend when present)
%
% Outputs
%   Mq : 7x7 mass matrix
%   Vv : 7x1 Coriolis/centrifugal vector, [] when dq is empty
%   Gv : 7x1 gravity vector
%
% See also CH3_CONTROL_AFFINE, CH3_IMPACT, CH3_ZD_POINT, CH3_MODEL_SIGNATURE.

lam = [];
if isfield(p, 'model_blend') && ~isempty(p.model_blend)
    lam = p.model_blend;
    if ~(isscalar(lam) && isreal(lam) && isfinite(lam) && lam >= 0 && lam <= 1)
        error('ch3_mvg:blend', ...
              'p.model_blend must be empty or a scalar in [0, 1] (got %s).', ...
              mat2str(lam));
    end
    if lam == 1, lam = []; end              % exactly today's model
end

need_V = nargout > 1 && ~isempty(dq);

if isempty(lam)
    Mq = M(q);
    if need_V, Vv = V([q; dq]); else, Vv = []; end
    if nargout > 2, Gv = G(q); end
else
    Mq = (1 - lam) * M_m30(q) + lam * M(q);
    if need_V
        x  = [q; dq];
        Vv = (1 - lam) * V_m30(x) + lam * V(x);
    else
        Vv = [];
    end
    if nargout > 2
        Gv = (1 - lam) * G_m30(q) + lam * G(q);
    end
end

end
