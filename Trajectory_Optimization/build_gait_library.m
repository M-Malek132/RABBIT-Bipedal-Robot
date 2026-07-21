%% =========================================================
% BUILD GAIT LIBRARY — speed continuation
%
% Marches the NEC1 desired-speed equality (p.v_des) outward from a seed
% gait's natural speed, in small steps, warm-starting each solve from the
% previous CONVERGED gait. Small steps + warm start are what keep the
% periodic orbit intact: jumping v_des far in one solve satisfies the speed
% but breaks periodicity (we saw 0.355 -> 0.5 blow periodicity up to 0.80).
%
% Produces a struct array `gait_library`, one entry per speed, and saves it.
% This is the raw material for gait_library.png / gait_interpolation.png.
%% =========================================================
clear; clc;

[p, lb, ub, options] = hzd_problem_setup();
options = optimoptions(options, 'Display', 'final');   % quieter per-solve

%% ---- Seed gait (a known-good periodic gait) ----
seed_file = 'hzd_result_2026-07-20_12-06-36.mat';   % natural speed ~0.355 m/s
seed_path = fullfile(fileparts(mfilename('fullpath')), '..', 'Results', seed_file);
S      = load(seed_path);
z_seed = S.z_opt;

% Measure the seed's natural speed so the continuation is centred on it.
p.v_des = 0.355;                       % placeholder so inspect can evaluate
rep0    = inspect_solution(z_seed, p);
v_nat   = rep0.walking_speed;
fprintf('\nSeed natural speed = %.4f m/s\n', v_nat);

%% ---- Speed schedule: centre on v_nat, march both ways ----
dv     = 0.02;                         % continuation step [m/s]
v_hi   = 0.45;                         % upper target
v_lo   = 0.25;                         % lower target
v_up   = (v_nat + dv) : dv : v_hi;     % speeds above natural
v_down = (v_nat - dv) :-dv : v_lo;     % speeds below natural

feas_tol = options.ConstraintTolerance;

gait_library = struct('v_des',{},'v_actual',{},'z_opt',{},'fval',{}, ...
                      'exitflag',{},'feasible',{},'max_periodicity',{}, ...
                      'step_length',{},'step_duration',{},'cost_norm',{},'report',{});

% Seed entry (already converged at its natural speed)
gait_library = add_entry(gait_library, v_nat, z_seed, S.fval, 1, rep0, feas_tol);

%% ---- Continue UP then DOWN, each warm-started from the previous gait ----
gait_library = continuation(gait_library, v_up,   z_seed, p, lb, ub, options, feas_tol, 'UP');
gait_library = continuation(gait_library, v_down, z_seed, p, lb, ub, options, feas_tol, 'DOWN');

%% ---- Sort by speed and report ----
[~, order]   = sort([gait_library.v_des]);
gait_library = gait_library(order);

fprintf('\n===================== GAIT LIBRARY =====================\n');
fprintf('%8s %10s %8s %10s %10s %12s %6s\n', ...
        'v_des','v_actual','L_step','T_step','maxPeriod','cost/L','feas');
for k = 1:numel(gait_library)
    g = gait_library(k);
    fprintf('%8.3f %10.4f %8.3f %10.4f %10.4f %12.1f %6s\n', ...
            g.v_des, g.v_actual, g.step_length, g.step_duration, ...
            g.max_periodicity, g.cost_norm, tf(g.feasible));
end
fprintf('=======================================================\n');

%% ---- Save the library ----
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
ts = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
lib_file = fullfile(results_dir, ['gait_library_', ts, '.mat']);
save(lib_file, 'gait_library', 'p');
fprintf('Gait library saved to: %s\n', lib_file);

%% =========================================================
% helpers
%% =========================================================
function lib = continuation(lib, v_schedule, z_start, p, lb, ub, options, feas_tol, tag)
    z_prev = z_start;
    for v = v_schedule
        p.v_des = v;
        [z_opt, fval, exitflag] = fmincon(@(z) hzd_cost(z, p), z_prev, ...
            [], [], [], [], lb, ub, @(z) hzd_constraints(z, p), options);

        rep = inspect_solution(z_opt, p);
        feasible = solution_is_feasible(rep, feas_tol);
        fprintf('[%s] v_des=%.3f -> v=%.4f  maxPeriod=%.4f  feasible=%s\n', ...
                tag, v, rep.walking_speed, max(abs(rep.periodicity)), tf(feasible));

        lib = add_entry(lib, v, z_opt, fval, exitflag, rep, feas_tol);

        if ~feasible
            fprintf('[%s] stopping: gait no longer feasible at v_des=%.3f\n', tag, v);
            break;   % do NOT warm-start from an infeasible gait
        end
        z_prev = z_opt;   % warm-start next speed from this converged gait
    end
end

function lib = add_entry(lib, v_des, z_opt, fval, exitflag, rep, feas_tol)
    lib(end+1) = struct( ...
        'v_des',           v_des, ...
        'v_actual',        rep.walking_speed, ...
        'z_opt',           z_opt, ...
        'fval',            fval, ...
        'exitflag',        exitflag, ...
        'feasible',        solution_is_feasible(rep, feas_tol), ...
        'max_periodicity', max(abs(rep.periodicity)), ...
        'step_length',     rep.step_length, ...
        'step_duration',   rep.step_duration, ...
        'cost_norm',       rep.cost_normalized, ...
        'report',          rep);
end

function ok = solution_is_feasible(rep, feas_tol)
    % Feasible = sim succeeded AND every constraint within tolerance.
    viol = max([ abs(rep.periodicity(:)); abs(rep.foot_height); ...
                 abs(rep.foot_velocity(:)); abs(rep.theta_end); ...
                 abs(rep.velocity_ceq); ...
                 max(0, rep.ground_c); max(0, rep.clearance_c); ...
                 max(0, rep.step_length_c) ]);
    ok = (rep.status > 0) && (viol <= feas_tol);
end

function s = tf(b)
    if b, s = 'YES'; else, s = 'no'; end
end
