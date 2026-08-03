function mu = ch5_min_norm_mu(psi, LgV, ny)
%CH5_MIN_NORM_MU  The Chapter-3 min-norm CLF control, in closed form.
%
%   mu = ch5_min_norm_mu(psi, LgV, ny)
%
%       mu = 0                          if psi <= 0
%       mu = -(psi / ||LgV||^2) LgV'    otherwise
%
% the least-norm point on the halfspace  LgV mu <= -psi.
%
% Used in two places, for two different reasons:
%
%   1. as the ANSWER, in ch5_ctrl_clf_qp with no input box;
%   2. as the WARM START for every constrained QP in this chapter.
%
% The second is not an optimization. quadprog started from the origin was
% failing outright -- exitflag -8 ("no step direction"), -3, or hitting the
% iteration cap -- on ~1% of pendulum samples, at states where the barrier row
% was perfectly well conditioned (|L_b| ~ 2, h ~ 1). The cause is distance: on
% this maneuver psi reaches ~1e5, so the CLF row alone demands ||mu|| of order
% 1e3, and an interior-point method asked to walk that far from z = 0 while a
% 1e6 slack penalty deforms the objective does not always get there.
%
% Starting from this point instead, the solver only has to correct for the
% barrier row -- which is usually inactive, and when active moves mu by a
% modest amount. The failures go away, and where they never occurred the
% answer is unchanged: it is the same QP, solved from a better place.
%
% See also CH5_CTRL_CLF_QP, CH5_CTRL_ECBF_CLF_QP, CH5_CTRL_CBF_CLF_QP.

nrm2 = LgV * LgV.';

if psi <= 0 || nrm2 < 1e-14
    mu = zeros(ny, 1);
else
    mu = -(psi / nrm2) * LgV.';
end

end
