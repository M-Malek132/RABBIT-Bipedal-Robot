function T = ch5_report(runs, varargin)
%CH5_REPORT  Tabulate a set of Chapter-5 runs against the chapter's claims.
%
%   T = ch5_report(runs)
%   T = ch5_report(runs, 'title', 'Serial spring-mass')
%
% runs is a struct array with fields .name and .sim (a ch5_simulate output).
%
% ------------------------------------------------------------- what is reported
% The column that matters is MIN h. Chapter 5's claim is forward invariance of
% C = {h >= 0}, so a single negative value falsifies it for that run, and no
% amount of good behaviour elsewhere compensates. It is therefore reported as a
% number rather than a tick, because HOW negative separates two very different
% situations:
%
%   * order 1e-1  the constraint was simply not enforced -- the baseline.
%   * order dt    the constraint WAS enforced, at sample instants, and the
%                 plant moved between them. Sampled-data enforcement of a
%                 continuous-time condition cannot do better than O(dt), and
%                 ch5_system's control_dt notes measure that scaling directly.
%
% MIN y_rb is the same claim one level up, per Remark 5.6: it is the value of
% the outermost set function y_rb of (5.28), which is what the ECBF row
% actually constrains. If min h is fine but min y_rb is not, the barrier row
% was violated somewhere and the constraint survived by luck.
%
% It is reported RELATIVE to Kb*eta_b, and that is not a cosmetic choice. The
% row is y_rb = mu_b + Kb eta_b >= 0, and on the pendulum Kb eta_b reaches
% several thousand; a raw y_rb of -3e-5 there is the solver resting on the
% constraint to its own tolerance, which is what an active constraint is
% supposed to look like. Judging it against zero in absolute terms would report
% a violation on every run in which the barrier ever did anything.
%
% %ACT is how much of the run the safety constraint was binding. It reads as
% the price of safety: on the pendulum at p2min = -0.5 it is roughly two-thirds
% of the run, which is the quantitative form of Fig. 5.4c's "the controller has
% to aggressively move the links".
%
% %INF counts samples where the QP had no solution. For an ECBF run with no
% input box this can only happen where the barrier row loses its grip on mu
% (L_g L_f^(rb-1) h * Ainv = 0), which on the pendulum is the arm fully
% extended. It is a genuine gap in the guarantee and is counted rather than
% smoothed over.
%
% Output
%   T : the same numbers as a struct array, for programmatic use
%
% See also CH5_SIMULATE, CH5_MAIN, CH5_PLOT_SPRINGMASS, CH5_PLOT_PENDULUM.

o = struct('title', '');
for k = 1:2:numel(varargin), o.(lower(varargin{k})) = varargin{k+1}; end

if ~isempty(o.title)
    fprintf('\n%s\n', o.title);
    fprintf('%s\n', repmat('=', 1, max(60, numel(o.title))));
end

fprintf('%-22s %10s %11s %11s %8s %8s %6s %6s %6s\n', ...
        'run', 'min h', 'safe?', 'min y_rb/s', 'max|u|', 'p99|u|', ...
        '%act', '%inf', 'cpu');
fprintf('%s\n', repmat('-', 1, 96));

T = struct('name', {}, 'h_min', {}, 'safe', {}, 'y_rb_rel', {}, ...
           'u_max', {}, 'u_p99', {}, 'pct_active', {}, ...
           'pct_infeasible', {}, 'cpu', {});

for i = 1:numel(runs)
    s = runs(i).sim;

    ok  = ~isnan(s.y_rb);
    if any(ok)
        scale = max(1, abs(s.Kb_eta(ok)));
        yr_rel = min(s.y_rb(ok) ./ scale);
    else
        yr_rel = NaN;
    end

    au = max(abs(s.u), [], 1);

    rec = struct( ...
        'name',           runs(i).name, ...
        'h_min',          s.h_min, ...
        'safe',           s.h_min >= 0, ...
        'y_rb_rel',       yr_rel, ...
        'u_max',          max(au), ...
        'u_p99',          ch5_prctile(au, 99), ...
        'pct_active',     100*mean(s.cbf_active), ...
        'pct_infeasible', 100*mean(~s.feasible), ...
        'cpu',            s.cpu);

    fprintf('%-22s %+10.2e %11s %11s %8.2f %8.2f %5.1f%% %5.2f%% %5.0fs\n', ...
            rec.name, rec.h_min, safe_tag(s), fmt(yr_rel), rec.u_max, ...
            rec.u_p99, rec.pct_active, rec.pct_infeasible, rec.cpu);

    % max|u| next to p99|u| because the min-norm CLF-QP applies
    % ||mu|| = psi/||LgV||, which spikes wherever LgV passes near zero. On the
    % pendulum baseline that is one sample in 30001 -- 435 Nm against a p99 of
    % 27 -- and quoting only the max would describe the controller by its
    % single worst instant.
    if rec.u_max > 5 * max(rec.u_p99, eps)
        [~, k] = max(au);
        fprintf('%24s(max|u| is a %d-sample spike at t = %.3f s; ', ...
                '', sum(au > 5*rec.u_p99), s.t(k));
        fprintf('min-norm mu = psi/||LgV||)\n');
    end

    if ~s.ok
        fprintf('%26s(run stopped: %s)\n', '', s.reason);
    end
    if ~isempty(s.adm) && ~s.adm.ok
        fprintf('%26s(Cor 5.2: %s)\n', '', s.adm.msg);
    end

    T(end+1) = rec; %#ok<AGROW>
end

fprintf('%s\n', repmat('-', 1, 96));

end

% ---------------------------------------------------------------------------
function s = safe_tag(sim)
%SAFE_TAG  Distinguish "not enforced" from "enforced, sampled".
%
% The threshold is the control period, not an arbitrary epsilon: an excursion
% that scales with dt is the discretization, and one that does not is the
% controller.
if sim.h_min >= 0
    s = 'SAFE';
elseif abs(sim.h_min) <= 10 * sim.p.control_dt
    s = 'safe~O(dt)';
else
    s = 'VIOLATED';
end
end

function s = fmt(v)
if isnan(v), s = '-'; else, s = sprintf('%+.3e', v); end
end
