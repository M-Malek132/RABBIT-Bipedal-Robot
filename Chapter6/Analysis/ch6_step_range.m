function R = ch6_step_range(x0, alpha, p, l_d)
%CH6_STEP_RANGE  Measure the step-length range one nominal gait can deliver.
%
%   R = ch6_step_range(x0, alpha, p)
%   R = ch6_step_range(x0, alpha, p, l_d)
%
% Section 6.1.4's headline number: "we achieve the step length variation of
% [25:60] cm" from a single nominal gait of 45 cm, and Section 6.4 exists
% because that range is limited. This measures the equivalent number for THIS
% robot and THIS gait instead of quoting the thesis's.
%
% For each desired foothold l_d it runs ONE step from the nominal fixed point
% onto a stone of width p.mc.stone_sz centred at l_d, and records three things
% that are usually conflated:
%
%   .landed     the foot came down inside the window
%   .feasible   AND the QP was feasible at every control sample
%   .clean      AND no enabled physical limit was violated
%
% They are separated because they fail in that order and mean different things.
% "landed but infeasible" is the interesting case: the placement was achieved,
% but during some samples the barrier row and the torque box had no common
% solution, so the ECBF guarantee lapsed and the controller fell back to
% minimising the violation. The foot arrived; the certificate did not.
%
% ---------------------------------------------------- one step, not thirty
% A single step from the periodic orbit is the BEST case: the state is exactly
% on the nominal fixed point, so nothing but the stone is asking the controller
% to deviate. A thirty-step run over random stones (which is what Fig. 6.5
% shows) starts each step from wherever the previous one ended, so its range is
% narrower. This measurement is therefore an upper bound on what the multi-step
% runs can do, and reading it as the multi-step answer would flatter the method.
%
% Inputs
%   x0, alpha : the nominal gait's fixed point and coefficients
%   p         : parameter struct
%   l_d       : desired footholds to try (default 0.15 : 0.025 : 0.60)
%
% Output
%   R : struct .l_d .l_s .landed .feasible .clean .u_max .n_infeasible
%              .range_landed .range_clean  (as [lo hi], NaN if empty)
%              .nominal      the unconstrained step length, for the % figures
%
% See also CH6_STEP, CH6_SIMULATE, CH6_TABLE61.

if nargin < 4 || isempty(l_d), l_d = 0.15 : 0.025 : 0.60; end

n = numel(l_d);
R = struct('l_d', l_d, 'l_s', nan(1,n), ...
           'landed', false(1,n), 'feasible', false(1,n), 'clean', false(1,n), ...
           'u_max', nan(1,n), 'n_infeasible', nan(1,n));

% The unconstrained step, for the percentages.
q0 = p;  q0.controller = 'iolin_pd';  q0.cbf.problem = 'none';
s0 = ch6_step(x0, alpha, q0);
R.nominal = s0.l_s;

fprintf('  nominal (no barrier) step length %.4f m; stone size %.0f cm\n', ...
        R.nominal, 100*p.mc.stone_sz);
fprintf('  %-7s %-9s %-8s %-9s %-8s %s\n', ...
        'l_d', 'l_s', 'landed', 'feasible', 'max|u|', 'infeas');

for i = 1:n
    q = p;
    q.cbf.problem = 'stones';
    q.stone = ch6_resolve_stone(p.stones, ...
                                l_d(i) - p.mc.stone_sz/2, ...
                                l_d(i) + p.mc.stone_sz/2);
    try
        s = ch6_step(x0, alpha, q);
    catch
        continue;
    end
    if ~s.ok, continue; end

    R.l_s(i)          = s.l_s;
    R.landed(i)       = s.in_window;
    R.n_infeasible(i) = sum(~s.log.feasible);
    R.feasible(i)     = s.in_window && R.n_infeasible(i) == 0;
    R.u_max(i)        = max(abs(s.log.u(:)));
    R.clean(i)        = R.feasible(i) && limits_ok(s, p);

    fprintf('  %-7.3f %-9.4f %-8s %-9s %-8.1f %d\n', ...
            l_d(i), s.l_s, yn(R.landed(i)), yn(R.feasible(i)), ...
            R.u_max(i), R.n_infeasible(i));
end

R.range_landed = span(l_d, R.landed);
R.range_clean  = span(l_d, R.clean);

fprintf('\n  landed in window over %s\n', fmt(R.range_landed, R.nominal));
fprintf('  and clean (feasible + limits) over %s\n', fmt(R.range_clean, R.nominal));

end

% ---------------------------------------------------------------------------
function ok = limits_ok(s, p)
ok = true;
if p.limits.enable.grf,    ok = ok && min(s.log.Fz) >= p.limits.Fz_min - 1e-6; end
if p.limits.enable.torque, ok = ok && max(abs(s.log.u(:))) <= p.limits.u_max + 1e-6; end
if p.limits.enable.friction
    m = s.log.mu_fric(s.log.Fz > 1e-6);
    ok = ok && (isempty(m) || max(m) <= p.limits.mu_s + 1e-6);
end
end

function r = span(l_d, mask)
i = find(mask);
if isempty(i), r = [NaN NaN]; else, r = [l_d(i(1)), l_d(i(end))]; end
end

function s = fmt(r, nom)
if any(isnan(r))
    s = 'NOTHING';
else
    s = sprintf('[%.0f : %.0f] cm  (%+.0f%% to %+.0f%% of the %.0f cm nominal)', ...
                100*r(1), 100*r(2), 100*(r(1)/nom - 1), 100*(r(2)/nom - 1), ...
                100*nom);
end
end

function s = yn(b)
if b, s = 'yes'; else, s = 'no'; end
end
