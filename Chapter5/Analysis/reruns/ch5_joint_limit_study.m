function ch5_joint_limit_study()
%CH5_JOINT_LIMIT_STUDY  Write the self-collision down, and see if both barriers hold.
%
%   ch5_joint_limit_study()
%
% The chapter's pendulum result enforces the end-effector height perfectly by
% folding the arm to theta2 = -3.10 rad (level -0.5 m) or -2.23 rad (-1.0 m):
% almost completely back on itself, a self-collision on hardware that the
% model does not contain. The chapter states the lesson -- a constraint that is
% not written down is not enforced -- and this demonstrates it, and meets the
% question Chapter 6 will face on the walking robot: can SEVERAL
% relative-degree-4 barriers be satisfied at the same time?
%
% Every run is the chapter's pendulum mission (theta1 from -pi to pi,
% theta2 -> 0, 15 s) under the ECBF-CLF-QP with the height barrier, and in
% addition a joint-limit barrier theta2 >= limit (p.ecbf.extra, same poles as
% the height row). Grid: level -1.0 and -0.5 m x limit none / -2.5 / -2.2 rad,
% on both solvers the chapter uses. At theta1 = 0 the arm must fold to
% |theta2| >= acos(-1 - level) to clear the height (1.57 rad at -1.0 m, 2.09
% at -0.5 m), so both limits leave a feasible band: [-2.5, -2.09] and
% [-2.2, -1.57] where they bind.
%
% Reported per run: completed or not, min height margin h, min theta2 and the
% joint barrier's margin, the share of samples each barrier row is active and
% BOTH are, QP feasibility, peak torque and |theta1(T) - pi|.
%
% Output: Results/reruns/ch5_joint_limit/ (runs.mat, joint_limit.log).
% Runtime: 12 runs of 20-60 s, about 10 minutes.
%
% See also CH5_EXTRA_BARRIERS, CH5_CTRL_ECBF_CLF_QP, CH5_BARRIER.

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
OUT  = fullfile(ROOT, 'Results', 'reruns', 'ch5_joint_limit');
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'joint_limit.log');
if exist(logf, 'file'), delete(logf); end

levels  = [-1.0 -0.5];
limits  = [NaN -2.5 -2.2];
solvers = {'rk4', 'ode45'};
R = struct('level', {}, 'limit', {}, 'solver', {}, 'ok', {}, 'reason', {}, ...
           'h_min', {}, 'th2_min', {}, 'hj_min', {}, 'act_h', {}, 'act_j', {}, ...
           'act_both', {}, 'feasible', {}, 'tau_max', {}, 'th1_err', {}, 'cpu', {});
logline(logf, '=== ch5_joint_limit_study | %s', datestr(now)); %#ok<TNOW1,DATST>

for lv = levels
    for jl = limits
        for s = solvers
            args = {'system', 'pendulum', 'controller', 'ecbfclfqp', ...
                    'constraint', lv, 'integrator', s{1}};
            if ~isnan(jl)
                args = [args, {'ecbf.extra', struct('type', 'theta_min', 'joint', 2, ...
                                                     'value', jl, 'poles', [])}]; %#ok<AGROW>
            end
            p = ch5_params(args{:});
            t0 = tic;
            r = struct('level', lv, 'limit', jl, 'solver', s{1}, 'ok', false, ...
                       'reason', '', 'h_min', NaN, 'th2_min', NaN, 'hj_min', NaN, ...
                       'act_h', NaN, 'act_j', NaN, 'act_both', NaN, 'feasible', false, ...
                       'tau_max', NaN, 'th1_err', NaN, 'cpu', NaN);
            try
                sim = ch5_simulate(p);
                n = sim.n;
                r.ok = sim.ok; r.reason = sim.reason;
                r.h_min = sim.h_min;
                r.th2_min = min(sim.x(2, 1:n));
                r.act_h = 100 * mean(sim.cbf_active(1:n));
                if ~isempty(sim.h_extra)
                    r.hj_min = sim.h_extra_min(1);
                    aj = sim.act_extra(1, 1:n);
                    r.act_j = 100 * mean(aj);
                    r.act_both = 100 * mean(aj & sim.cbf_active(1:n));
                end
                r.feasible = all(sim.feasible(1:n));
                r.tau_max = max(abs(sim.u(:, 1:n)), [], 'all');
                r.th1_err = abs(sim.x(1, n) - pi);
            catch err
                r.reason = ['error: ' err.message];
            end
            r.cpu = toc(t0);
            R(end+1) = r; %#ok<AGROW>
            logline(logf, ['level %+.1f | limit %5s | %-5s | ok %d | h_min %+.4f | ' ...
                           'min theta2 %+.3f | joint margin %+.4f | active h %4.1f%% ' ...
                           'j %4.1f%% both %4.1f%% | QP feasible %d | max|tau| %7.1f | ' ...
                           '|theta1(T)-pi| %.4f | %s [%.0f s]'], lv, num2str(jl), s{1}, ...
                    r.ok, r.h_min, r.th2_min, r.hj_min, r.act_h, r.act_j, r.act_both, ...
                    r.feasible, r.tau_max, r.th1_err, r.reason, r.cpu);
        end
    end
end
save(fullfile(OUT, 'runs.mat'), 'R');
logline(logf, '=== JOINT_LIMIT_DONE');
fprintf('JOINT_LIMIT_DONE\n');
end

function logline(logf, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(logf, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
