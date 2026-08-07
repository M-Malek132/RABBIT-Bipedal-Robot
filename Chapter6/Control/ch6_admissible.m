function a = ch6_admissible(x, alpha, p, t)
%CH6_ADMISSIBLE  Is this step's initial condition inside the barrier's sets?
%
%   a = ch6_admissible(x, alpha, p, t)
%
% Corollary 5.2: the ECBF poles are admissible for a given x0 only if
%
%       p_i >= -ydot_(i-1)(x0) / y_(i-1)(x0),
%
% which for the rb = 2 chain of Chapter 6 is the pair of conditions
%
%       g(x0)     >= 0            the state starts in the safe set
%       h_CBF(x0) >= 0            i.e.  gdot(x0) >= -gamma_b g(x0)
%
% The barrier guarantees forward invariance FROM the set, not convergence TO it.
% A step that begins with h_CBF < 0 has no guarantee at all, and the run will
% then look like a controller failure when it is an initialisation failure.
% Checking it takes one evaluation, so there is no reason to discover it from a
% violated constraint 200 ms later.
%
% ============================================== WHERE THE SWING FOOT ACTUALLY STARTS
% Worth being concrete, because the sign is counterintuitive and it decides
% whether g_ST2 starts positive or negative.
%
% At the start of step k the stance foot is the foot that just LANDED and the
% swing foot is the one it landed past, so
%
%       l_f(0) = -L_(k-1)          (BEHIND the stance foot)
%       l_f(T) = +L_k              = the step length at impact       (6.13)
%
% and the swing foot travels through the stance foot, not away from it.
% Measured on the reference gait: l_f(0) = -0.353 m, g_ST1 = 0.953,
% g_ST2 = 0.315 -- both comfortably positive, because at l_f < 0 the foot is
% outside circle O2 on its LEFT branch.
%
% This is why "outside O2" is a real constraint rather than a formality. O2's
% left edge sits at l_f = l_min/2 - rho2, a few centimetres behind the stance
% foot, so the foot must climb OVER the arch to reach the l_f >= l_min branch
% where it is allowed to land. It cannot slide along the ground from one branch
% to the other. That climb is the scuffing guarantee of Remark 6.2, and it is
% enforced exactly over the final approach -- the part of the swing where
% scuffing actually happens.
%
% ------------------------------------------------------------ what can still fail
% Admissibility is not automatic. g_ST1 = R1 + l_max - |r - (-R1,0)| shrinks as
% the foot moves BACKWARD too, so a gait whose previous step was longer than
% l_max can start the next step already outside circle O1 -- the barrier would
% then be asked to recover rather than to maintain, which is not what it
% guarantees. That is a genuine limitation of using one nominal gait over a
% wide range of stone spacings, and it is one of the things the gait library of
% Section 6.4 fixes.
%
% Inputs
%   x     : 14x1 start-of-step state
%   alpha : ny x n_ctrl coefficients (unused, kept for signature symmetry with
%           the controllers -- barriers do not depend on the gait)
%   p     : parameter struct
%   t     : time offset to evaluate the stone at (default 0)
%
% Output
%   a : struct
%       .ok      all barriers admissible
%       .g       1 x nb   barrier values at x0
%       .h_cbf   1 x nb   gamma_b g + gdot at x0
%       .margin  1 x nb   min(g, h_cbf), the distance to inadmissibility
%       .labels  1 x nb cell
%       .why     text description when ~ok
%
% See also CH6_BARRIER, CH6_CBF_ROW, CH5_ECBF_ADMISSIBLE, CH6_SIMULATE.

if nargin < 4 || isempty(t), t = 0; end
if nargin < 2, alpha = []; end %#ok<NASGU>

[B, ~] = ch6_barrier(x, [], p, t);
nb = numel(B);

a = struct('ok', true, 'g', nan(1,nb), 'h_cbf', nan(1,nb), ...
           'margin', nan(1,nb), 'labels', {cell(1,nb)}, 'why', '');

bad = {};
for i = 1:nb
    g  = B(i).g;
    hc = p.cbf.gamma_b * g + B(i).gdot;

    a.g(i)      = g;
    a.h_cbf(i)  = hc;
    a.margin(i) = min(g, hc);
    a.labels{i} = B(i).label;

    if a.margin(i) < 0
        a.ok = false;
        bad{end+1} = sprintf('%s (g = %.3e, h_CBF = %.3e)', ...
                             B(i).label, g, hc);   %#ok<AGROW>
    end
end

if ~a.ok
    a.why = ['initial condition outside the barrier sets: ' ...
             strjoin(bad, '; ')];
end

end
