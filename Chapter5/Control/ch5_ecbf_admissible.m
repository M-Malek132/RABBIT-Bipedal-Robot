function adm = ch5_ecbf_admissible(x0, p, e)
%CH5_ECBF_ADMISSIBLE  Corollary 5.2's initial-condition test on the poles.
%
%   adm = ch5_ecbf_admissible(x0, p)
%   adm = ch5_ecbf_admissible(x0, p, e)      with a prebuilt ch5_ecbf_gain
%
% Theorem 5.1 does NOT say "pick stable poles and you are safe". It says C_0 is
% forward invariant whenever x0 lies in EVERY set of the chain,
%
%       C_i = { x : y_i(x) >= 0 },   y_i = (d/dt + p_i) o ... o (d/dt + p_1) h
%
% for i = 0, ..., rb. Only y_0 = h >= 0 is the constraint the user asked for;
% the other rb conditions are extra, and they are conditions on x0 AND the
% poles together. Corollary 5.2 turns them into the pole rule
%
%       p_i >= -ydot_(i-1)(x0) / y_(i-1)(x0),   i = 1, ..., rb.       (Cor 5.2)
%
% -------------------------------------------------------- what this means in use
% Read the rule for what it is: if the barrier is already heading toward the
% boundary at x0, the poles must be FAST ENOUGH to have arrested it in time.
% Slow poles are not a conservative choice here -- they are an inadmissible one.
% This is the limitation Section 5.2.3 closes with ("the choice of pole location
% depends on initial conditions ... requiring careful choice of these poles"),
% and it is checkable, so it is checked, at the actual x0 the simulation starts
% from rather than in a footnote.
%
% ------------------------------------------------ exactness, stated honestly
% y_i for i <= rb-1 uses derivatives h^(0..rb-1), all of which are in eta_b and
% all of which are exact functions of x0.
%
% The LAST rung is different. y_rb and ydot_(rb-1) need h^(rb), and h^(rb)
% contains the control -- which at t = 0 is the output of the very QP whose
% feasibility is in question. So those two are evaluated at the DRIFT value
% L_f^rb h, i.e. at u = 0, and reported separately as .top_is_drift. The drift
% value is the right diagnostic anyway: y_rb(x0) >= 0 under drift means the
% barrier row is already satisfied with no control at all, and y_rb(x0) < 0
% means the QP must be active from the first sample.
%
% Inputs
%   x0 : nx x 1 initial state
%   p  : parameter struct
%   e  : optional, from ch5_ecbf_gain (rebuilt if omitted)
%
% Output
%   adm : struct with
%           .ok            all exactly-checkable conditions hold
%           .y             (rb+1) x 1, y_0(x0) ... y_rb(x0)
%           .y_ok          (rb+1) x 1 logical, y_i >= 0
%           .p_min         rb x 1, the Corollary 5.2 lower bound per pole
%           .poles .slack  poles - p_min
%           .binding       index of the tightest condition
%           .top_is_drift  true (see above)
%           .msg           a one-line human summary
%
% See also CH5_ECBF_GAIN, CH5_BARRIER, CH5_SIMULATE.

b  = ch5_barrier(x0, p);
rb = b.rb;

if nargin < 3 || isempty(e)
    e = ch5_ecbf_gain(p, rb);
end
poles = e.poles(:);

% hd(m) = h^(m-1).  The top entry is the drift value, see the header.
hd = [b.eta_b(:); b.Lfrb];

%% ---------------------------------- the family (5.28), built by its own recursion
% y_i = sum_j c_i(j+1) h^(i-j) with c_i = poly(-p_1..p_i). Building it from the
% polynomial rather than by iterating (5.30) means an error in one is not
% mirrored in the other -- ch5_test_ecbf checks the two constructions agree.
y = zeros(rb+1, 1);
for i = 0:rb
    c = poly(-poles(1:i));                      % [1, a1, ..., ai]
    acc = 0;
    for j = 0:i
        acc = acc + c(j+1) * hd(i-j+1);
    end
    y(i+1) = acc;
end

%% -------------------------------------------------- Corollary 5.2, pole by pole
% ydot_(i-1) is y_(i-1) with every derivative index shifted up one.
p_min = -inf(rb, 1);
for i = 1:rb
    c = poly(-poles(1:i-1));
    ydot_prev = 0;
    for j = 0:(i-1)
        ydot_prev = ydot_prev + c(j+1) * hd(i-j+1);
    end
    y_prev = y(i);

    if y_prev > 0
        p_min(i) = -ydot_prev / y_prev;
    elseif y_prev == 0
        % On the boundary of C_(i-1): no finite pole rescues a strictly
        % inward-pointing violation, and any pole works if ydot >= 0.
        if ydot_prev >= 0, p_min(i) = 0; else, p_min(i) = inf; end
    else
        p_min(i) = inf;                          % x0 is already outside C_(i-1)
    end
end

slack = poles - p_min;

% Only i = 0..rb-1 are exactly checkable; the top rung is drift-evaluated.
y_ok = y >= 0;
ok   = all(y_ok(1:rb)) && all(slack(1:rb-1) >= -1e-9);

[~, binding] = min(slack);

if ok
    msg = sprintf(['poles admissible at x0 (tightest: p_%d = %.3f vs required ' ...
                   '%.3f)'], binding, poles(binding), p_min(binding));
else
    bad = find(~y_ok(1:rb), 1) - 1;
    if isempty(bad)
        msg = sprintf('pole p_%d = %.3f is below the required %.3f', ...
                      binding, poles(binding), p_min(binding));
    else
        msg = sprintf(['x0 is outside C_%d (y_%d = %.4f < 0): Theorem 5.1 ' ...
                       'needs x0 in every C_i, so forward invariance of ' ...
                       'C_0 does not follow'], bad, bad, y(bad+1));
    end
end

adm = struct('ok', ok, 'y', y, 'y_ok', y_ok, 'p_min', p_min, ...
             'poles', poles, 'slack', slack, 'binding', binding, ...
             'top_is_drift', true, 'rb', rb, 'msg', msg);

end
