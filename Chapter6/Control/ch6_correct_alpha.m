function [alpha, info] = ch6_correct_alpha(x, alpha, p)
%CH6_CORRECT_ALPHA  Start the virtual constraints where the robot actually is.
%
%   [alpha_c, info] = ch6_correct_alpha(x_plus, alpha, p)
%
% Rewrites the first two Bezier columns so that at the start-of-step state
% x_plus the outputs and their rates are ZERO:
%
%       y(x_plus) = 0,   ydot(x_plus) = 0
%
% and leaves every other column alone. This is the post-impact correction of
% Westervelt et al., "Feedback Control of Dynamic Bipedal Robot Locomotion",
% Sec. 6.6, used whenever a step ends somewhere the gait did not plan for.
%
% WHY CHAPTER 6 NEEDS IT. A gait's virtual constraints are designed so that
% Delta maps its own pre-impact state onto y = 0, ydot = 0. A step the barrier
% has LENGTHENED OR SHORTENED is not that state: MEASURED on posture_195, a step
% that lands at 0.461 m instead of 0.427 m leaves |y| = 0.23 rad and |ydot| =
% 4.4 rad/s after impact. Driving that out needs more than the robot has --
% unconstrained PD asks for 1138 Nm and a NEGATIVE normal force (-124 N) -- so
% under the 250 Nm box and the contact rows the QP cannot, the torso pitches
% forward, and the swing foot comes down 6 cm in front of the stance foot. That
% happens with the barrier rows removed, so it is not a barrier failure; it is
% what walking a fixed gait off its orbit costs. Every multi-step stepping-stone
% run failed at step 2 for this reason, even with every stone within 2.5 cm of
% the nominal step.
%
% WHAT IT CHANGES. For a Bezier polynomial of degree M in s,
%
%       b(s) = sum_k alpha_k B_k(s),
%
% b and b' at a phase s0 are linear in (alpha_0, alpha_1), so matching
% b(s0) = H q+ and b'(s0) = H dq+ / sdot+ is a 2x2 solve per output. At
% s0 = 0 it is the textbook pair alpha_0 = H q+, alpha_1 = alpha_0 +
% H dq+ / (M sdot+). The correction's weight, B_0 + B_1, falls below 10% by
% s = 0.4 for M = 5, so the second half of the step is the gait as designed and
% the impact at the end of it is the gait's own.
%
% It is applied to EVERY step by ch6_simulate (p.post_impact_correction), for
% every controller Table 6.1 compares, and on the periodic orbit it is a no-op
% to the residual of the orbit itself.
%
% It is skipped, and alpha returned unchanged, when the correction is not well
% defined: a basis other than Bezier, a phase outside [-0.25, 0.25] of the
% step start, or a phase rate too small to divide by.
%
% Inputs
%   x     : 14x1 start-of-step state
%   alpha : ny x (M+1) Bezier coefficients
%   p     : parameter struct (gait fields from ch6_load_gait)
%
% Outputs
%   alpha : corrected coefficients
%   info  : .applied .s0 .sdot .y0 .ydot0 (before) .why (when skipped)
%
% See also CH6_SIMULATE, CH3_OUTPUTS, CH3_BEZIER.

[y0, yd0, o] = ch3_outputs(x, alpha, p);
info = struct('applied', false, 's0', o.s, 'sdot', o.sdot, ...
              'y0', y0, 'ydot0', yd0, 'why', '');

if ~strcmpi(p.basis, 'bezier')
    info.why = 'basis is not bezier';            return;
end
if abs(o.s) > 0.25
    info.why = sprintf('phase %.3f too far from the step start', o.s);  return;
end
if abs(o.sdot) < 0.5
    info.why = sprintf('phase rate %.3f too small', o.sdot);           return;
end

M  = size(alpha, 2) - 1;
s0 = o.s;
q  = x(1:p.nq);
dq = x(p.nq+1:2*p.nq);

% Bernstein values and slopes of B_0, B_1 at s0.
B  = [(1-s0)^M,          M*s0*(1-s0)^(M-1)];
dB = [-M*(1-s0)^(M-1),   M*(1-s0)^(M-1) - M*(M-1)*s0*(1-s0)^(M-2)];
A  = [B; dB];

% What the untouched columns already contribute at s0.
rest = alpha;  rest(:, 1:2) = 0;
[b_rest, db_rest] = ch3_bezier(rest, s0);

target   = p.H * q;                  % b(s0)
target_d = (p.H * dq) / o.sdot;      % b'(s0) = d(Hq)/ds

for i = 1:size(alpha, 1)
    alpha(i, 1:2) = (A \ [target(i) - b_rest(i); target_d(i) - db_rest(i)]).';
end
info.applied = true;
end
