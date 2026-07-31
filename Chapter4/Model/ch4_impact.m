function [x_plus, impulse] = ch4_impact(x_minus, p, unc)
%CH4_IMPACT  The reset map Delta for the TRUE (perturbed) model.
%
%   [x_plus, impulse] = ch4_impact(x_minus, p)        uses p.uncertainty
%   [x_plus, impulse] = ch4_impact(x_minus, p, unc)   given perturbation
%   [x_plus, impulse] = ch4_impact(x_minus, p, [])    nominal -> ch3_impact
%
% Same three steps as ch3_impact -- rigid plastic impact, relabel, re-plant --
% but with the perturbed mass matrix in the momentum balance.
%
% THE TWO PERTURBATIONS BEHAVE COMPLETELY DIFFERENTLY HERE, and it is worth
% being precise about why.
%
%   mass_scale.  The balance is
%
%       [ sM   -Jsw' ] [ dq+     ]   [ sM dq- ]
%       [ Jsw    0   ] [ impulse ] = [   0    ]
%
%   Substituting impulse = s*w makes s divide out of the first row entirely and
%   the second row never had it, so dq+ is EXACTLY INDEPENDENT of s and only
%   the impulse scales, linearly. A uniformly heavier robot lands with the same
%   post-impact velocity and a proportionally larger contact impulse. This is
%   checked numerically in ch4_test_model, not just asserted -- and it is the
%   reason Case IV can scale the robot by 3 without the impact map exploding.
%
%   load_mass.  A point mass at the torso base changes M NON-UNIFORMLY (only
%   the two base entries), so nothing divides out: the post-impact velocity
%   genuinely changes. Physically the extra torso inertia resists the impulse
%   that the leg collision tries to transmit through the body.
%
% Because dq+ is scale-invariant, calling ch3_impact for a pure mass_scale
% would give the right STATE and the wrong IMPULSE. The impulse is a reported
% Table 3.1 quantity, so this function computes it properly rather than
% reusing a map that happens to agree on the part we look at most.
%
% Inputs
%   x_minus : 14x1 pre-impact state
%   p       : parameter struct
%   unc     : perturbation; omitted -> p.uncertainty; [] -> nominal
%
% Outputs
%   x_plus  : 14x1 start-of-next-step state
%   impulse : 2x1 contact impulse [Ns] for the TRUE model
%
% See also CH3_IMPACT, CH4_CONTROL_AFFINE, CH3_RELABEL.

if nargin < 3, unc = p.uncertainty; end

s  = 1; mL = 0;
if ~isempty(unc)
    if isfield(unc, 'mass_scale') && ~isempty(unc.mass_scale), s  = unc.mass_scale; end
    if isfield(unc, 'load_mass')  && ~isempty(unc.load_mass),  mL = unc.load_mass;  end
end

if s == 1 && mL == 0
    [x_plus, impulse] = ch3_impact(x_minus, p);
    return;
end

nq = p.nq;
q  = x_minus(1:nq);
dq = x_minus(nq+1:2*nq);

% --- 1. plastic impact on the perturbed inertia --------------------------
M_mat = s * M(q);
if mL ~= 0
    M_mat(1,1) = M_mat(1,1) + mL;
    M_mat(2,2) = M_mat(2,2) + mL;
end

Jsw = J_sw(q);                   % PRE-impact swing-foot Jacobian
nc  = size(Jsw, 1);

A = [M_mat, -Jsw.'; ...
     Jsw,    zeros(nc)];
b = [M_mat * dq; ...
     zeros(nc, 1)];

sol     = A \ b;
dq_plus = sol(1:nq);
impulse = sol(nq+1:end);

x_post = [q; dq_plus];           % q is continuous through impact

% --- 2. relabel -----------------------------------------------------------
x_rel = ch3_relabel(x_post, p);

% --- 3. re-plant the new stance foot on the ground ------------------------
foot   = P_st(x_rel(1:nq));
x_plus = x_rel;
x_plus(2) = x_plus(2) + foot(2);

end
