function ch3_speed_wall()
%CH3_SPEED_WALL  Which limit stops the speed march at 1.15 m/s?
%
% ch3_speed_ladder marches 1.5628 -> 1.15 in 0.05 steps and then refuses every
% further step: 0.05, 0.025 and 0.013 all return in ~15 s with the NEC1
% equality untouched, so it is not a step-size problem. At 1.15 three rows are
% active simultaneously -- peak torque on its 195 Nm box, stance friction at
% exactly mu_s = 0.400, and the normal force exactly on its 50 N floor -- and
% no feasible direction is left that also lowers the speed.
%
% This asks 1.10 m/s four times, each with ONE limit relaxed, and reports which
% relaxation lets the march continue. A run that takes says that limit is the
% wall; a run that refuses in ~15 s says it is not. Nothing here is a proposal
% to change the limits -- it is a diagnostic to find which one is binding.
%
% Output: Results/reruns/speed_ladder/wall_probe.log

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'speed_ladder');
logf = fullfile(SPD, 'wall_probe.log');

L  = load(fullfile(SPD, 'gait_v1150.mat'));
p0 = L.p; z0 = L.z;
p0.v_des = 1.10;

probes = {
  'baseline (nothing relaxed)', @(q) q
  'torque box 195 -> 230 Nm',   @(q) setfield(q, 'limits', setfield(q.limits, 'u_max',   230))
  'friction mu_s 0.40 -> 0.60', @(q) setfield(q, 'limits', setfield(q.limits, 'mu_s',    0.60))
  'Fz floor 50 -> 10 N',        @(q) setfield(q, 'limits', setfield(q.limits, 'Fz_min',  10))
};

logln(logf, '=== wall probe at v_des 1.10 from the 1.15 gait | %s', datestr(now));
for i = 1:size(probes, 1)
    [name, f] = probes{i, :};
    p = f(p0);
    t0 = tic;
    [z, out] = ch3_col_solve(p, z0, struct('MaxFunctionEvaluations', 3.0e5, ...
                                           'MaxIterations', 400));
    secs = toc(t0);
    E = ch3_col_eval(z, p);
    took = out.max_ceq <= 1e-3;
    logln(logf, '%-28s -> %s | speed %.4f | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
          name, ternary(took, 'TOOK   ', 'refused'), E.L_step / E.T, ...
          out.max_c, out.max_ceq, secs);
    if took
        save(fullfile(SPD, sprintf('wall_probe_%d.mat', i)), 'z', 'p', 'out');
    end
end
logln(logf, '=== WALL_PROBE_DONE');
fprintf('WALL_PROBE_DONE\n');
end

function s = ternary(b, x, y)
if b, s = x; else, s = y; end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
