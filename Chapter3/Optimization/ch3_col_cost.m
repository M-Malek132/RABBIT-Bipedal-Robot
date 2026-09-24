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

% ------------------------------------------------ energy per STEP, not per metre
% p.cost_normalize = false drops the division and minimizes int ||u||^2 alone.
%
% WHAT IT ACTUALLY CHANGES. Unnormalized, a shorter step is simply cheaper, so
% the stride goes to whatever lower bound holds it -- the row 3 floor, 0.15 m
% by default. That is not the step-in-place degeneracy the header warns about,
% which the floor already excludes; it is a gait pinned to its shortest legal
% stride. Useful when the stride is being MARCHED, where the row 19 ceiling
% sets the target and the floor is far below it, and wrong otherwise.
%
% WHY THE MARCH NEEDS IT. Normalized, every millimetre of stride removed is a
% cost increase, so a solve asked to shorten the stride has no descent
% direction and simply refuses.
%
% MEASURED. Asked to bring the stride from 0.4376 m under a 0.40 m ceiling, SQP
% stalled at a step length of 1e-7 with feasibility frozen at 0.0376, exactly
% the violation. Under a 0.43 m ceiling -- a violation of 0.0076 -- it stalled
% in two iterations with the stride unmoved. The fight is structural, not a
% matter of rung size.
%
% DEFAULT TRUE, so every gait solved so far means what it meant. Turn it off
% only with the ceiling on, or the step-in-place degeneracy comes back.
if isfield(p, 'cost_normalize') && ~p.cost_normalize
    if ~(isfield(p.limits, 'enable') && isfield(p.limits.enable, 'step_len_max') ...
            && p.limits.enable.step_len_max)
        error('ch3_col_cost:unboxedStride', ...
              ['p.cost_normalize = false without the step-length ceiling ' ...
               '(p.limits.enable.step_len_max): int ||u||^2 alone is minimized by ' ...
               'the shortest legal stride, so the gait collapses onto the row 3 ' ...
               'floor (%.2f m) rather than onto a stride anyone chose. Enable the ' ...
               'ceiling and march it, or leave the objective normalized.'], ...
              p.step_len_min);
    end
    J = cost_scale(p) * E.int_u2;
    return;
end

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

J = cost_scale(p) * E.int_u2 / L_eff;

end

% ---------------------------------------------------------------------------
function s = cost_scale(p)
%COST_SCALE  Multiplier on the objective. Default 1; the optimum is unchanged
% by any positive value, but the MERIT function's balance is not.
%
% WHY IT IS NEEDED. fmincon weighs the objective against the constraint
% violation. Here the objective is ~6.8e3 with a gradient norm of ~3.6e4, while
% the step-length row's gradient is order 1, so a violated stride ceiling is
% nearly invisible next to the cost. MEASURED on posture_195, asked to bring
% the stride in by 1 mm: it moved 1.3e-5 m in 83 iterations and 30 minutes,
% creeping in the right direction while the solver spent itself on the cost.
% Asked for 0.1 mm it stalled outright. Turning ScaleProblem off instead let it
% move but diverged to a spurious solution 36 off a true rollout.
%
% Shrinking the objective makes feasibility the cheaper thing to buy. Use it
% for a march whose constraint the objective opposes; leave it at 1 otherwise,
% since a tiny objective also means a loose optimum.
s = 1;
if isfield(p, 'cost_scale') && ~isempty(p.cost_scale)
    s = p.cost_scale;
    if ~(isscalar(s) && isfinite(s) && s > 0)
        error('ch3_col_cost:scale', ...
              'p.cost_scale must be a positive scalar (got %s).', mat2str(s));
    end
end
end
