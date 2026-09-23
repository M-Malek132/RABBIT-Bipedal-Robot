function ch3_torque_wall()
%CH3_TORQUE_WALL  What stops the torque march at 162.5 Nm?
%
% ch3_torque_march took the box 195 -> 180 -> 170 -> 165 -> 162.5, every rung
% landing and the last one cleanly (max|c| 7.3e-11, verifies at 7.1e-05). Then
% 160 refused twice, at residuals of 6.2e-04 and 9.2e-04 on the torque row --
% four orders worse than the rung above, which says barrier rather than step
% size. At every rung two other rows are exactly active: stance friction at
% mu_s = 0.400 and the normal force on its 50 N floor.
%
% This asks for 160 Nm from the 162.5 gait five times, each with ONE thing
% changed, and reports which change lets the box come down. It is a diagnostic,
% not a proposal to move the limits: the point is to name what the declared
% 120 Nm box would cost, since ch3_params sets that limit and Chapter 4's
% sec:box120 shows every applied-torque law dying under it.
%
% The last probe is the interesting one. Peak torque fell as the gait got
% FASTER along this march (1.2815 -> 1.2921 m/s), so pinning the speed DOWN may
% make the box harder rather than easier -- worth knowing before anyone tries
% to hold 1.2 m/s and 120 Nm at once.
%
% Output: Results/reruns/torque_march/wall_probe.log

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'torque_march');
logf = fullfile(SPD, 'wall_probe.log');

L  = load(fullfile(SPD, 'gait_u162.mat'));
z0 = L.z_try;
p0 = L.p;
p0.limits.u_max = 160;

probes = {
  'nothing relaxed',            @(q) q
  'friction mu_s 0.40 -> 0.50', @(q) setfield(q, 'limits', setfield(q.limits, 'mu_s',   0.50))
  'Fz floor 50 -> 25 N',        @(q) setfield(q, 'limits', setfield(q.limits, 'Fz_min', 25))
  'both of the above',          @(q) setfield(q, 'limits', setfield(setfield(q.limits, 'mu_s', 0.50), 'Fz_min', 25))
  'speed pinned at 1.20 m/s',   @(q) pin(q, 1.20)
};

logln(logf, '=== torque wall probe at u_max 160 from the 162.5 gait | %s', datestr(now));
for i = 1:size(probes, 1)
    [name, f] = probes{i, :};
    p = f(p0);
    t0 = tic;
    [z, out] = ch3_col_solve(p, z0, struct('MaxFunctionEvaluations', 4e5, ...
                                           'MaxIterations', 400));
    E = ch3_col_eval(z, p);
    V = ch3_col_verify(z, p, false);
    lam = [E.lam, E.lamm];
    took = out.max_c <= 1e-4 && out.max_ceq <= 1e-6 && V.ok;
    logln(logf, ['%-28s -> %s | peak|u| %6.1f | v %.4f | mu %.3f | minFz %5.1f | ' ...
                 'max|c| %.2e | max|ceq| %.2e | verify %d | %.0f s'], ...
          name, ternary(took, 'TOOK   ', 'refused'), ...
          max([abs(E.u(:)); abs(E.um(:))]), E.L_step/E.T, ...
          max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), ...
          out.max_c, out.max_ceq, V.ok, toc(t0));
    if took
        save(fullfile(SPD, sprintf('wall_probe_%d.mat', i)), 'z', 'p', 'out');
    end
end
logln(logf, '=== TORQUE_WALL_DONE');
fprintf('TORQUE_WALL_DONE\n');
end

function q = pin(q, v)
q.enforce_nec1 = true;
q.v_des = v;
end

function s = ternary(b, x, y)
if b, s = x; else, s = y; end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
