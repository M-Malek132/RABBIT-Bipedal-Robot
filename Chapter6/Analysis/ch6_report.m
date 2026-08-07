function R = ch6_report(sim, p)
%CH6_REPORT  What a stepping-stone run actually did, printed and returned.
%
%   R = ch6_report(sim, p)
%
% Follows the Chapter-3 rule: EVERY quantity is printed with its measured value
% whether or not it is enforced, and [E] marks the enforced ones. Measuring
% before tightening is the only workflow that produces limits a solve can meet,
% and it is the reason ch3_report exists in the form it does.
%
% What gets printed, and why each line is there:
%
%   FOOTSTEP PLACEMENT   the headline of Figs. 6.5, 6.10, 6.25a -- step length
%     at impact against the window that was live AT IMPACT. For a moving stone
%     that is not the window the step started with.
%
%   BARRIERS             min over the run of g_i and of the row margin. g_i < 0
%     means the constraint was violated; margin < 0 means the ECBF condition was
%     violated, which happens FIRST and is the earlier warning. Both are
%     reported because they fail in that order.
%
%   QP HEALTH            infeasible samples, and samples where a barrier row had
%     no grip on u. Remark 6.6 attributes part of Table 6.1's residual failure
%     rate to infeasibility, so it is counted rather than inferred.
%
%   PHYSICAL LIMITS      (6.25) and input saturation: normal force, friction
%     cone, torque. Measured always, enforced when p.limits.enable says so.
%
%   GAIT LIBRARY         how often (6.23) had to EXTRAPOLATE, which is Remark
%     6.6's other named cause.
%
% Inputs
%   sim : output of ch6_simulate
%   p   : parameter struct
%
% Output
%   R : struct with every printed number, so a sweep can collect them
%
% See also CH6_SIMULATE, CH6_PLOT_STONES, CH6_TABLE61, CH3_REPORT.

R = struct();
n = sim.n_ok;

fprintf('\n================= CHAPTER 6 REPORT =================\n');
fprintf(' controller %s | barrier %s (%s) | poles (%.3g, %.3g)\n', ...
        p.controller, p.cbf.problem, p.cbf.form, p.cbf.gamma_b, p.cbf.gamma);
fprintf(' %d/%d steps completed, stones move: %s\n', n, sim.n_req, p.stones.motion);

if n == 0
    fprintf(' nothing to report: %s\n', sim.reason);
    R.n_ok = 0; R.success = false; R.reason = sim.reason;
    fprintf('====================================================\n');
    return;
end

%% ------------------------------------------------------ footstep placement
fprintf('\n FOOTSTEP PLACEMENT  (l_min <= l_s <= l_max at impact)     (6.7)\n');
fprintf('   %-4s %-9s %-9s %-9s %-9s %s\n', ...
        'step', 'l_min', 'l_s', 'l_max', 'slack', 'in');
for k = 1:n
    slack = min(sim.l_s(k) - sim.l_min(k), sim.l_max(k) - sim.l_s(k));
    fprintf('   %-4d %-9.4f %-9.4f %-9.4f %-+9.4f %s\n', ...
            k, sim.l_min(k), sim.l_s(k), sim.l_max(k), slack, ...
            yn(sim.hit(k)));
end
R.l_s   = sim.l_s;
R.l_min = sim.l_min;
R.l_max = sim.l_max;
R.hit   = sim.hit;
R.n_hit = sum(sim.hit);
fprintf('   %d/%d placements inside the window\n', R.n_hit, n);

%% ------------------------------------------------------------- the barriers
nb = 0;
for k = 1:n, nb = max(nb, size(sim.steps(k).log.h, 1)); end

if nb > 0
    fprintf('\n BARRIERS   g >= 0 is the constraint, margin >= 0 is the ECBF row\n');
    g_min = inf(1, nb);  m_min = inf(1, nb);  n_act = zeros(1, nb);
    for k = 1:n
        lg = sim.steps(k).log;
        g_min = min(g_min, min(lg.h, [], 2).');
        m_min = min(m_min, min(lg.margin, [], 2).');
        n_act = n_act + sum(lg.cbf_active, 2).';
    end
    B0 = ch6_barrier(sim.steps(1).x(:,1), [], setstone(p, sim.steps(1)), 0);
    for j = 1:nb
        fprintf('   %-22s min g %+9.4f   min margin %+10.3e   active %d samples\n', ...
                B0(j).label, g_min(j), m_min(j), n_act(j));
    end
    R.g_min = g_min;  R.margin_min = m_min;  R.n_active = n_act;
end

%% -------------------------------------------------------------- QP health
n_inf  = 0;  n_samp = 0;  n_dead = 0;  d_max = 0;
for k = 1:n
    lg = sim.steps(k).log;
    n_inf  = n_inf  + sum(~lg.feasible);
    n_samp = n_samp + numel(lg.feasible);
    n_dead = n_dead + sum(lg.n_dead > 0);
    d_max  = max(d_max, max(lg.delta));
end
fprintf('\n QP HEALTH\n');
fprintf('   infeasible samples        %d / %d\n', n_inf, n_samp);
fprintf('   rows with no grip on u    %d samples\n', n_dead);
fprintf('   max CLF slack delta       %.4g  (0 = the exponential bound held)\n', d_max);
R.n_infeasible = n_inf;  R.n_samples = n_samp;  R.n_dead = n_dead;
R.delta_max = d_max;

%% -------------------------------------------------------- physical limits
Fz = []; mu = []; uu = [];
for k = 1:n
    lg = sim.steps(k).log;
    Fz = [Fz, lg.Fz];                        %#ok<AGROW>
    mu = [mu, lg.mu_fric(lg.Fz > 1e-6)];     %#ok<AGROW>
    uu = [uu, lg.u];                         %#ok<AGROW>
end
fprintf('\n PHYSICAL LIMITS   [E] = enforced in the QP                (6.25)\n');
pr('normal force  Fz_min',  min(Fz),      p.limits.Fz_min, p.limits.enable.grf,      'ge');
pr('friction      |Fx/Fz|', max(mu),      p.limits.mu_s,   p.limits.enable.friction, 'le');
pr('torque        max |u|', max(abs(uu(:))), p.limits.u_max, p.limits.enable.torque, 'le');
R.Fz_min = min(Fz);  R.mu_max = max(mu);  R.u_max = max(abs(uu(:)));

%% ---------------------------------------------------------- gait library
if any(sim.lib_used)
    fprintf('\n GAIT LIBRARY                                        (6.22)-(6.24)\n');
    fprintf('   steps taking an interpolated gait   %d / %d\n', sum(sim.lib_used), n);
    fprintf('   steps needing EXTRAPOLATION         %d   (Remark 6.6)\n', ...
            sum(sim.lib_extrap));
    R.n_extrap = sum(sim.lib_extrap);
end

%% ----------------------------------------------------------------- verdict
fprintf('\n VERDICT  %s   (%s)\n', upper(yn(sim.success)), sim.reason);
fprintf('====================================================\n');

R.n_ok    = n;
R.success = sim.success;
R.reason  = sim.reason;

end

% ---------------------------------------------------------------------------
function pr(name, val, lim, enabled, sense)
switch sense
    case 'ge', ok = val >= lim - 1e-6;
    case 'le', ok = val <= lim + 1e-6;
end
tag = '   ';
if enabled, tag = '[E]'; end
fprintf('   %s %-22s %10.4g   limit %-10.4g %s\n', tag, name, val, lim, ...
        pf(ok, enabled));
end

function s = pf(ok, enabled)
if ~enabled
    if ok, s = '(ok, not enforced)'; else, s = '(EXCEEDED, not enforced)'; end
else
    if ok, s = 'ok'; else, s = '*** VIOLATED ***'; end
end
end

function s = yn(b)
if b, s = 'yes'; else, s = 'no'; end
end

function q = setstone(p, step)
q = p;
q.stone = step.stone;
end
