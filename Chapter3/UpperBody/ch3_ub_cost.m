function J = ch3_ub_cost(p)
%CH3_UB_COST  Effort penalty on an OPTIMIZED upper-body wrench.
%
%   J = ch3_ub_cost(p)
%
%       J = w * mean_s || W(s) ./ scale ||^2 ,   s on a 41-point grid in [0,1]
%
% Zero unless p.ub.optimize is set. Without this term a free external wrench
% costs nothing, so the optimizer would lean on it to erase the joint torque
% integral entirely; the weight w is the exchange rate between "Newtons the
% support has to provide" and "joint torque the legs have to provide", and it
% is the knob a crutch/assistance study sweeps.
%
% The phase average is state-independent (W is a function of s only), which is
% why this can be added to the stage-3 objective without touching ch3_col_eval.
%
% See also CH3_UB_DEFAULTS, CH3_COL_COST.

[~, optimize] = ch3_ub_info(p);
J = 0;
if ~optimize || p.ub.w == 0
    return;
end
sg  = linspace(0, 1, 41);
acc = 0;
for k = 1:numel(sg)
    Wk  = ch3_bezier(p.ub.beta, sg(k)) ./ p.ub.scale;
    acc = acc + sum(Wk.^2);
end
J = p.ub.w * acc / numel(sg);
end
