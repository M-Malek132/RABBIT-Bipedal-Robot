function p = ch3_col_effective_params(z, p)
%CH3_COL_EFFECTIVE_PARAMS  p with the phase endpoints the decision vector holds.
%
%   p = ch3_col_effective_params(z, p)
%
% With p.free_theta false this returns p unchanged. With it true, theta_minus
% and theta_plus live in z (ch3_col_pack), and every function that evaluates
% the gait -- the phase s(q), the outputs, the zero dynamics, the node-1 and
% node-N equalities, the verification rollout -- must read THOSE, not the
% values p happened to carry when the solve started. This is the one place
% that copies them across, so the collocation, the check and the saved gait
% cannot disagree about what s = 0 and s = 1 mean.
%
% A gait solved with free theta is saved with this p (ch3_col_solve returns it
% as out.p), so a later consumer that reads p.theta_minus / p.theta_plus --
% ch4_load_gait, the controllers, ch3_simulate -- gets the solved values.
%
% See also CH3_COL_UNPACK, CH3_COL_EVAL, CH3_COL_SOLVE.

if ~(isfield(p, 'free_theta') && ~isempty(p.free_theta) && p.free_theta)
    return;
end

[~, ~, ~, theta_pm] = ch3_col_unpack(z, p);
if ~(theta_pm(2) > theta_pm(1))
    error('ch3_col_effective_params:order', ...
          ['theta_plus (%.4f) must exceed theta_minus (%.4f): s = (theta - ' ...
           'theta_minus)/(theta_plus - theta_minus) is undefined otherwise. ' ...
           'Keep p.theta_bounds disjoint.'], theta_pm(2), theta_pm(1));
end
p.theta_minus = theta_pm(1);
p.theta_plus  = theta_pm(2);

end
