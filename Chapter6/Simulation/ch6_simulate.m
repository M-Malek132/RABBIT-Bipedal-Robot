function out = ch6_simulate(x0, gait, p, terrain)
%CH6_SIMULATE  Walk a set of discrete footholds.  Figs. 6.5, 6.6, 6.10, 6.25.
%
%   out = ch6_simulate(x0, gait, p)              nominal window every step
%   out = ch6_simulate(x0, gait, p, terrain)     a terrain from ch6_terrain
%
% One step per stone. Before each step the stone is resolved into p.stone (the
% radii become numbers), the gait for the step is chosen, and ch6_step runs the
% closed loop to strike. After the step the placement is scored against the
% window AS IT STOOD AT IMPACT.
%
% ------------------------------------------------------------------ gait choice
% `gait` is either
%
%   a numeric ny x n_ctrl matrix -- ONE nominal gait for every step. This is
%     Sections 6.1-6.3 and controller II of (6.27); the chapter is explicit that
%     "for the work of using CBF for stepping stones, we use only one nominal
%     walking gait", and that this is what limits the step-length range.
%
%   a gait library struct from ch6_lib_load -- Section 6.4. alpha is
%     interpolated per step from the stone's desired step length, (6.22)-(6.23),
%     and the CBF then handles the transient the interpolation does not.
%
% Both go through the same code path so Table 6.1's columns differ only in what
% they are supposed to differ in.
%
% ------------------------------------------------------ what counts as a failure
% A trial fails on ANY of:
%
%   * the guard never fires (the robot did not complete the step)
%   * the state goes non-finite
%   * the step length at impact lands outside the foothold window
%   * an ENABLED physical limit was violated (friction, normal force, torque)
%
% -- which is the thesis's own criterion: "the controller is considered
% successful for a trial run if the bipedal robot is able to walk over this
% terrain without violation of foot placement, ground contact, friction, and
% input constraints".
%
% An INFEASIBLE QP is deliberately not on that list. It is the mechanism by
% which those get violated, not a fifth criterion, and since the controller
% keeps the contact rows and the input box hard even when it cannot meet the
% barriers, an infeasible sample frequently violates nothing. Counting it would
% penalise a controller for reporting honestly that its guarantee lapsed while
% rewarding one that silently saturated. ch6_report counts and prints it
% separately; Remark 6.6 is about that number.
%
% Inputs
%   x0      : 14x1 start state (post-impact, stance foot on the ground)
%   gait    : ny x n_ctrl alpha, or a library struct
%   p       : parameter struct
%   terrain : optional 1 x n struct array from ch6_terrain. Missing or empty
%             uses [p.stones.l_min, p.stones.l_max] for p.n_steps steps.
%
% Output
%   out : struct
%     .steps    1 x k struct array of ch6_step outputs, each with .stone and
%               .alpha added
%     .l_s      1 x k step length at impact [m]
%     .l_min .l_max   1 x k window at impact [m]
%     .hit      1 x k logical, placement inside the window
%     .n_ok     completed steps
%     .success  every attempted step completed AND satisfied everything
%     .failed   ~success
%     .reason   why it stopped or failed
%     .t .x     concatenated trajectory
%     .lib_used 1 x k logical, gait came from the library
%     .lib_extrap 1 x k logical, the interpolation extrapolated       (Rmk 6.6)
%
% See also CH6_STEP, CH6_TERRAIN, CH6_LIB_ALPHA, CH6_TABLE61, CH6_REPORT.

if nargin < 4 || isempty(terrain)
    terrain = repmat(struct('l_d',   0.5*(p.stones.l_min + p.stones.l_max), ...
                            'l_min', p.stones.l_min, ...
                            'l_max', p.stones.l_max), 1, p.n_steps);
end

n = numel(terrain);
use_lib = isstruct(gait) && isfield(gait, 'alpha');

steps = {};
t_all = []; x_all = []; t_off = 0;
x     = x0(:);

l_s   = nan(1, n);  l_lo = nan(1, n);  l_hi = nan(1, n);
hit   = false(1, n);
extrap = false(1, n);

success = true;
reason  = 'completed';

for k = 1:n

    if ~all(isfinite(x))
        success = false;
        reason  = sprintf('non-finite state entering step %d', k);
        break;
    end

    pk = p;
    pk.stone = ch6_resolve_stone(p.stones, terrain(k).l_min, terrain(k).l_max);

    if use_lib
        [alpha_k, li] = ch6_lib_alpha(gait, terrain(k).l_d);
        extrap(k) = li.extrapolated;
    else
        alpha_k = gait;
    end

    s = ch6_step(x, alpha_k, pk);
    s.stone = pk.stone;
    s.alpha = alpha_k;
    s.k     = k;

    steps{end+1} = s;                       %#ok<AGROW>
    t_all = [t_all, t_off + s.t];           %#ok<AGROW>
    x_all = [x_all, s.x];                   %#ok<AGROW>
    t_off = t_off + s.T;

    l_s(k)  = s.l_s;
    l_lo(k) = s.stone_end.l_min;
    l_hi(k) = s.stone_end.l_max;
    hit(k)  = s.in_window;

    if ~s.ok
        success = false;
        reason  = sprintf('step %d never reached the guard (T = %.3f s)', k, s.T);
        break;
    end

    v = check_step(s, pk);
    if ~isempty(v)
        success = false;
        reason  = sprintf('step %d: %s', k, v);
        break;
    end

    x = s.x_next;
end

if isempty(steps)
    steps_arr = struct([]);
else
    steps_arr = [steps{:}];
end

nk = numel(steps);

out = struct('steps', steps_arr, ...
             'l_s', l_s(1:nk), 'l_min', l_lo(1:nk), 'l_max', l_hi(1:nk), ...
             'hit', hit(1:nk), ...
             'n_ok', nk, 'n_req', n, ...
             'success', success && nk == n, 'failed', ~(success && nk == n), ...
             'reason', reason, ...
             't', t_all, 'x', x_all, 'x_final', x, ...
             'lib_used', repmat(use_lib, 1, nk), ...
             'lib_extrap', extrap(1:nk), ...
             'terrain', terrain);

end

% ---------------------------------------------------------------------------
function why = check_step(s, p)
%CHECK_STEP  The thesis's success criterion, applied to one completed step.
%
% Only ENABLED limits are checked. A disabled limit is still measured and
% reported by ch6_report -- the Chapter-3 measure-then-tighten workflow -- but
% it cannot fail a trial it was not asked to enforce, or Table 6.1 would be
% scoring the gait rather than the controller.
%
% AN INFEASIBLE QP IS NOT ITSELF A FAILURE HERE. The criterion is about the four
% quantities above, and the thesis says so: "without violation of foot
% placement, ground contact, friction, and input constraints". Infeasibility is
% the MECHANISM by which those get violated, not a fifth criterion -- and since
% ch6_ctrl_cbf_clf_qp keeps the contact rows and the box hard even when it
% cannot meet the barriers, an infeasible sample very often violates nothing at
% all. Counting it as a failure would score the controller down for reporting
% honestly that its guarantee lapsed, while a controller that silently saturated
% would score better. It is counted and printed by ch6_report instead.
why = '';

if ~s.in_window
    why = sprintf('foot placed at l_s = %.4f, outside [%.4f, %.4f]', ...
                  s.l_s, s.stone_end.l_min, s.stone_end.l_max);
    return;
end

if p.limits.enable.grf && any(s.log.Fz < p.limits.Fz_min - 1e-6)
    why = sprintf('normal force fell to %.1f N (floor %.1f N)', ...
                  min(s.log.Fz), p.limits.Fz_min);
    return;
end

if p.limits.enable.friction
    % Only meaningful where the foot is loaded; |Fx/Fz| at Fz ~ 0 is a division
    % by a vanishing number, not a slip. The GRF row is what keeps Fz away from
    % zero, and it is checked above.
    m = s.log.mu_fric(s.log.Fz > 1e-6);
    if any(m > p.limits.mu_s + 1e-6)
        why = sprintf('friction cone violated, |Fx/Fz| reached %.3f (limit %.3f)', ...
                      max(m), p.limits.mu_s);
        return;
    end
end

if p.limits.enable.torque && any(abs(s.log.u(:)) > p.limits.u_max + 1e-6)
    why = sprintf('torque reached %.1f Nm (limit %.1f Nm)', ...
                  max(abs(s.log.u(:))), p.limits.u_max);
    return;
end

end
