function ch3_mu_sens(mu_list, u_max)
%CH3_MU_SENS  What a tighter friction coefficient costs the reference gait.
%
%   ch3_mu_sens()                 the ladder 0.38 0.36 0.35
%   ch3_mu_sens([0.38 0.36 0.35])
%   ch3_mu_sens(0.35, 210)        one rung, with the torque box opened to 210
%
% posture_195 lands with its impact impulse EXACTLY on the mu = 0.4 cone
% (row 13 of ch3_col_constraints reads +5.6e-9). So the landing is only just
% consistent with the coefficient Table 3.1 assumes. This re-solves the gait
% at tighter mu and reports what the cone costs.
%
% A CONTINUATION LADDER, NOT A JUMP. Going straight from 0.40 to 0.35 stalls
% at exitflag -2 with the PEAK TORQUE row over by 1.3e-3 (and 2.4e-3 at 0.30):
% the cone itself ends up satisfied, but the landing it forces wants more than
% the 195 Nm box, which the baseline gait is already sitting exactly on. So
% each rung warm-starts from the rung before it rather than from the baseline,
% which is what unstuck the realizability march on this project.
%
% u_max, when given, opens the torque box for every rung -- the second way of
% answering the same question: let the box breathe and see whether the cone is
% merely expensive rather than infeasible.
%
% Scaling ON: with it off, SQP runs away from a gait it has already converged
% (measured on this project).
if nargin < 1 || isempty(mu_list), mu_list = [0.38 0.36 0.35]; end
if nargin < 2, u_max = []; end
ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'mu_sens');
if ~exist(SPD,'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'mu_sens.log');

S = load(fullfile(ROOT,'Results','ch3_gait_posture_195.mat'));
z0 = S.z; p0 = ch3_upgrade_params(S.p);
if ~isempty(u_max), p0.limits.u_max = u_max; end
tag = ''; if ~isempty(u_max), tag = sprintf('_u%03.0f', u_max); end
report(logf, p0, z0, p0.limits.mu_s, 'baseline', NaN);

z_prev = z0;                     % THE LADDER: each rung starts from the last
for mu = mu_list(:).'
    f = fullfile(SPD, sprintf('gait_mu%03.0f%s.mat', 100*mu, tag));
    if exist(f,'file')
        L = load(f); z_prev = L.z;
        report(logf, L.p, L.z, mu, 'resolved', L.out.wall_time);
        continue;
    end
    p = p0; p.limits.mu_s = mu; p.scale_problem = true;
    t0 = tic;
    [z, out] = ch3_col_solve(p, z_prev);
    out.wall_time = toc(t0);
    save(f, 'z', 'p', 'out');
    logln(logf, 'mu %.2f%s: exitflag %d, fval %.4f, max|c| %.2e, max|ceq| %.2e, %.0f s', ...
          mu, tag, out.exitflag, out.fval, out.max_c, out.max_ceq, out.wall_time);
    report(logf, p, z, mu, 'resolved', out.wall_time);
    if out.max_c > 1e-6
        logln(logf, '  rung mu %.2f is INFEASIBLE (max|c| %.2e) -- later rungs start from it anyway', ...
              mu, out.max_c);
    end
    z_prev = z;
end
fprintf('MU_SENS_DONE\n');
end

function report(logf, p, z, mu, tag, secs)
E = ch3_col_eval(z, p);
[c, ceq] = ch3_col_constraints(z, p);
V = ch3_col_verify(z, p, false);
pk = max([abs(E.u(:)); abs(E.um(:))]);
logln(logf, ['%-9s mu %.2f | peak|u| %7.2f Nm | T %.4f s | L %.4f m | v %.4f m/s | ' ...
             'cost %.4f | impulse-cone row %+.3e | NIC1 row %+.3e | NIC2 row %+.3e | ' ...
             'max|c| %+.2e | max|ceq| %.2e | verify %d (%.1e) | %.0f s'], ...
      tag, mu, pk, E.T, E.L_step, E.L_step/E.T, ch3_col_cost(z,p), ...
      c(13), c(6), c(5), max(c), max(abs(ceq)), V.ok, V.max_dev, secs);
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
