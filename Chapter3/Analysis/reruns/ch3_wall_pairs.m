function ch3_wall_pairs()
%CH3_WALL_PAIRS  Can any PAIR of relaxations get the march below 1.15 m/s?
%
% ch3_speed_wall showed no SINGLE relaxation releases the wall: at 1.15 the
% torque box, the stance friction cone and the Fz floor are all active at once,
% and relaxing any one of them still refuses in ~15 s. Three rows active means
% a feasible descent direction may need two of them to move together.
%
% Each probe asks for 1.10 m/s from the verified 1.15 gait. A probe that TAKES
% says that pair is the wall, and the march can continue with those two
% relaxed and re-tighten them later (the torque ladder in
% ch3_stage3_from_scratch is the pattern for re-tightening).
ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'speed_ladder');
logf = fullfile(SPD, 'wall_pairs.log');
L = load(fullfile(SPD, 'gait_v1150.mat'));
p0 = L.p; z0 = L.z; p0.v_des = 1.10;

set2 = @(q, a, va, b, vb) setfield(q, 'limits', ...
        setfield(setfield(q.limits, a, va), b, vb));
probes = {
  'torque 230 + Fz floor 10',   @(q) set2(q, 'u_max', 230, 'Fz_min', 10)
  'torque 230 + mu_s 0.60',     @(q) set2(q, 'u_max', 230, 'mu_s',   0.60)
  'mu_s 0.60 + Fz floor 10',    @(q) set2(q, 'mu_s',  0.60,'Fz_min', 10)
  'all three relaxed',          @(q) setfield(q, 'limits', setfield(set2(q, 'u_max', 230, 'Fz_min', 10).limits, 'mu_s', 0.60))
};
logln(logf, '=== wall pairs at v_des 1.10 from the 1.15 gait | %s', datestr(now));
for i = 1:size(probes,1)
    [name, f] = probes{i,:};
    p = f(p0);
    t0 = tic;
    [z, out] = ch3_col_solve(p, z0, struct('MaxFunctionEvaluations', 3e5, 'MaxIterations', 400));
    E = ch3_col_eval(z, p);
    took = out.max_ceq <= 1e-3;
    logln(logf, '%-26s -> %s | speed %.4f | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
          name, ternary(took,'TOOK   ','refused'), E.L_step/E.T, out.max_c, out.max_ceq, toc(t0));
    if took, save(fullfile(SPD, sprintf('pair_probe_%d.mat', i)), 'z', 'p', 'out'); end
end
logln(logf, '=== WALL_PAIRS_DONE');
fprintf('WALL_PAIRS_DONE\n');
end
function s = ternary(b,x,y), if b, s=x; else, s=y; end, end
function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:}); fprintf('%s\n', s);
fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
