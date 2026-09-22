function M = ch6_momentum(x0, alpha, p, l_d, poles)
%CH6_MOMENTUM  What a footstep barrier costs the NEXT step.
%
%   M = ch6_momentum(x0, alpha, p)
%   M = ch6_momentum(x0, alpha, p, l_d, poles)
%
% ch6_step_range answers "can one step land on the stone". On posture_195 the
% answer is yes over 0.20-0.575 m, and it is the wrong question, because a
% walking robot has to take the step after. This measures that.
%
% THE MECHANISM. The swing foot of posture_195 strikes at 7.5 m/s forward
% relative to the stance foot. The upper-bound row (6.9), h_CBF = gamma_b g +
% gdot >= 0, lets the foot close on the far edge of the stone no faster than
% gamma_b times the room left, so on a +-2.5 cm stone a pole of 30 starts
% braking the foot with about 25 cm to go -- long before the edge, even when
% the nominal step would have landed inside. The step stretches, the hip slows,
% and the underactuated robot arrives at the next step without the momentum to
% take it: the next step fails with the barrier removed and every limit in
% place, and succeeds only with them lifted (thousands of Nm, a pulling normal
% force). A pole high enough to leave the momentum alone acts only in the last
% centimetre and can barely move the foot.
%
% Two measurements:
%
%   .sweep    per desired foothold l_d (split poles from p): one step onto the
%             stone from the orbit, then one step with NO barrier from where it
%             left the robot (virtual constraints corrected, ch6_correct_alpha).
%             Records landed, step duration, hip speed at impact, and whether
%             the second step completed with its foot ahead of the stance foot.
%
%   .poles    the stone centred on the NOMINAL step, upper-bound pole swept
%             over `poles`: how much a barrier that has nothing to fix still
%             slows the robot.
%
% Inputs
%   x0, alpha : the nominal gait's fixed point and coefficients
%   p         : parameter struct (from ch6_load_gait)
%   l_d       : footholds for .sweep (default 0.20:0.05:0.60)
%   poles     : upper-bound poles for .poles (default [30 60 100 200 400])
%
% Output
%   M : struct
%     .nominal   .l_s .T .v_hip (no barrier)
%     .sweep     .l_d .l_s .landed .T .v_hip .next_ok .next_l_s
%     .poles     .pole .l_s .landed .T .v_hip .next_ok .next_l_s
%
% See also CH6_STEP_RANGE, CH6_CORRECT_ALPHA, CH6_SIMULATE.

if nargin < 4 || isempty(l_d),   l_d   = 0.20 : 0.05 : 0.60;     end
if nargin < 5 || isempty(poles), poles = [30 60 100 200 400];    end

VX = p.nq + 1;                          % horizontal velocity of the hip

q0 = p;  q0.cbf.problem = 'none';
s0 = ch6_step(x0, alpha, q0);
M.nominal = struct('l_s', s0.l_s, 'T', s0.T, 'v_hip', s0.x_end(VX));

sz = p.mc.stone_sz;
fprintf('  nominal step %.4f m in %.3f s, hip %.2f m/s at impact\n', ...
        s0.l_s, s0.T, s0.x_end(VX));

% ---------------------------------------------------------------- the sweep
n = numel(l_d);
S = struct('l_d', l_d, 'l_s', nan(1,n), 'landed', false(1,n), ...
           'T', nan(1,n), 'v_hip', nan(1,n), ...
           'next_ok', false(1,n), 'next_l_s', nan(1,n));
fprintf('  %-6s %-7s %-6s %-6s %-6s %s\n', 'l_d', 'l_s', 'in', 'T', 'v_hip', 'next step');
for i = 1:n
    q = p;
    q.stone = ch6_resolve_stone(p.stones, l_d(i) - sz/2, l_d(i) + sz/2);
    [S.l_s(i), S.landed(i), S.T(i), S.v_hip(i), S.next_ok(i), S.next_l_s(i)] = ...
        two_steps(x0, alpha, q, q0, VX);
    fprintf('  %-6.3f %-7.4f %-6d %-6.3f %-6.2f %d (%.3f)\n', l_d(i), S.l_s(i), ...
            S.landed(i), S.T(i), S.v_hip(i), S.next_ok(i), S.next_l_s(i));
end
M.sweep = S;

% ---------------------------------------------------------------- the poles
m = numel(poles);
P = struct('pole', poles, 'l_s', nan(1,m), 'landed', false(1,m), ...
           'T', nan(1,m), 'v_hip', nan(1,m), ...
           'next_ok', false(1,m), 'next_l_s', nan(1,m));
fprintf('  stone centred on the nominal step, upper-bound pole swept\n');
for i = 1:m
    q = p;
    q.cbf.poles_by_label = [p.cbf.poles_by_label; ...
                            {'ls <= lmax', poles(i) * [1 1]}];
    q.stone = ch6_resolve_stone(p.stones, s0.l_s - sz/2, s0.l_s + sz/2);
    [P.l_s(i), P.landed(i), P.T(i), P.v_hip(i), P.next_ok(i), P.next_l_s(i)] = ...
        two_steps(x0, alpha, q, q0, VX);
    fprintf('  pole %-4d l_s %.4f T %.3f v_hip %.2f next %d\n', poles(i), ...
            P.l_s(i), P.T(i), P.v_hip(i), P.next_ok(i));
end
M.poles = P;
end

% ---------------------------------------------------------------------------
function [l_s, landed, T, v, next_ok, next_l_s] = two_steps(x0, alpha, q, q0, VX)
l_s = NaN; landed = false; T = NaN; v = NaN; next_ok = false; next_l_s = NaN;
s1 = ch6_step(x0, alpha, q);
if ~s1.ok, return; end
l_s = s1.l_s;  landed = s1.in_window;  T = s1.T;  v = s1.x_end(VX);
if ~all(isfinite(s1.x_next)), return; end

ac = ch6_correct_alpha(s1.x_next, alpha, q0);
s2 = ch6_step(s1.x_next, ac, q0);
next_l_s = s2.l_s;
% "Completed" = the guard fired with the foot ahead of the stance foot by at
% least the chapter's smallest stone, and no enabled contact limit broke.
next_ok = s2.ok && s2.l_s > 0.15 && min(s2.log.Fz) >= q0.limits.Fz_min - 1e-6;
end
