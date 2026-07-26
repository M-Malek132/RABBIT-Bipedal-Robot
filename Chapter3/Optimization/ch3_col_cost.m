function J = ch3_col_cost(z, p)
%CH3_COL_COST  Stage 3 objective: torque-squared per unit distance.
%
%   J = ch3_col_cost(z, p)
%
% Chapter 3 minimizes
%
%       J = (1 / L_step) * int_0^T ||u(t)||_2^2 dt
%
% WHY THE DIVISION BY L_STEP IS NOT COSMETIC.  Without it, a shorter step is
% simply cheaper -- the integral shrinks with the distance travelled -- so the
% optimizer is actively rewarded for shrinking the stride until the robot
% steps in place. Normalizing turns the objective into ENERGY PER DISTANCE,
% which is the quantity anyone actually cares about, and removes the
% degeneracy at the source rather than patching it with a constraint.
%
% A step-length floor is still enforced as an inequality, but it is a backstop:
% the normalization is what makes the objective well posed.
%
% GUARDING THE DIVISION.  L_step can pass through zero while the optimizer is
% still far from feasible, which would make J infinite and hand fmincon an
% undefined gradient. So the denominator is floored SMOOTHLY: below the floor
% the cost transitions to a quadratic penalty rather than a cliff, keeping J
% continuous and differentiable everywhere. A discontinuous objective is far
% more damaging to an SQP method than a slightly inaccurate one.
%
% See also CH3_COL_EVAL, CH3_COL_CONSTRAINTS.

E = ch3_col_eval(z, p);

L_min = p.step_len_min;
L     = E.L_step;

if L >= L_min
    L_eff = L;
else
    % Smooth continuation below the floor: value and slope match at L = L_min,
    % and L_eff stays strictly positive for any L (including negative L, i.e.
    % a backwards step), so J never becomes Inf or NaN.
    L_eff = L_min / (1 + ((L_min - L)/L_min)^2);
end

J = E.int_u2 / L_eff;

end
