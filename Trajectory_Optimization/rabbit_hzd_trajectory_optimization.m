%% =========================================================
% PROBLEM SETUP — parameters, bounds, fmincon options
% (shared with build_gait_library.m via hzd_problem_setup.m)
%% =========================================================
[p, lb, ub, options] = hzd_problem_setup();

% Pin the walking speed for THIS single solve. The seed gait's measured
% natural speed is ~0.355 m/s, so v_des = 0.355 makes it an already-feasible
% fixed point (clean validation). To reach other speeds use build_gait_library.
p.v_des = 0.355;                % [m/s]

% ---- PHYSICAL CONSTRAINTS: enable ONE AT A TIME, re-solving each ----------
% Order matters. The measured reference gait has the vertical GRF going
% NEGATIVE (min ~ -193 N: the foot would have to be pulled down, i.e. it lifts
% off), which ALSO makes the friction ratio |Fx/Fz| blow up at the Fz=0
% crossing. So the ground-reaction constraint is the ROOT and must come first;
% friction only becomes meaningful once Fz stays positive. Impulse (~10 Ns) and
% mean GRF (~309 N) are already fine. Warm-start each phase from the previous:
%   phase 1: min vertical GRF   (keep the foot in compression, Fz > 0)  <== NOW
%   phase 2: + friction cone    (now well-defined)
%   phase 3: + torque           (seed peaks ~257 Nm > u_max 120)
%   phase 4: + impact impulse   (already ~10 Ns; likely free)
%   (then p.enforce_stability for the final, hardest phase)
p.enforce_grf      = true;      % <-- phase 1 (the root pathology)
p.enforce_friction = false;
p.enforce_torque   = false;
p.enforce_impulse  = false;

%% =========================================================
% INITIAL STATE — natural post-impact start-of-step pose
% (px pz qt q1 q2 q3 q4  dpx dpz dqt dq1 dq2 dq3 dq4)
% Constructed so that: stance foot on ground; SWING FOOT BEHIND the
% stance foot (-0.28 m) with 6 cm clearance; theta = -0.15 rad (hip
% behind stance foot, so phase sweeps -0.15 -> theta_plus=0.30 forward);
% both knees bent naturally; stance-foot velocity = 0 at t=0.
%% =========================================================
x0 = [ +0.0000; -0.9310; +0.2999; -0.7934; +0.6869; -0.6318; +0.9810; ...
       +0.3952; -0.0419; +0.1847; +0.1847; +0.1045; +0.0000; +0.0000];

r  = check_ground_validity(x0(1:7), x0(8:14));
if ~r.is_valid
    disp(r.violations);
end

%% =========================================================
% INITIAL SPLINE COEFFICIENTS  (see make_coeffs.m for rationale)
%% =========================================================
coeffs0 = make_coeffs(x0, p);

%% =========================================================
% SETTLE THE INITIAL STATE (Poincare / fixed-point iteration)
% DISABLED (n_settle = 0) -- tried and it makes things WORSE.
%
% Iterating the step map only converges if the gait is already near a
% STABLE periodic orbit. coeffs0 is a crude linear interpolation, far
% from periodic, so the iteration diverges instead of settling: theta
% blew up to +/-2..3 rad (robot flipping over) within 2 steps and the
% resulting x0 violated the bounds and ground constraints.
%
% Worth revisiting ONLY once fmincon can produce a near-periodic gait --
% then settling could polish it. Set n_settle > 0 to re-enable.
%% =========================================================
n_settle = 0;
if n_settle > 0
    [x0, settle_hist] = settle_initial_state(coeffs0, x0, p, n_settle);
    coeffs0 = make_coeffs(x0, p);   % rebuild so spline(s=0) == settled joints
end

z0 = [coeffs0(:); x0];

%% =========================================================
% WARM START (optional but strongly recommended once a gait exists)
% Set warm_start_file to a previously saved result to start fmincon from
% that gait instead of the cold x0. Cold starts from x0 are a gamble (even
% the first success wandered for 100+ iterations); warm starting from a
% known periodic gait is what makes speed continuation / a gait library
% work: solve one gait, nudge p.v_des, re-solve from the previous solution.
%
% The saved p may predate newer fields (e.g. v_des); we keep the CURRENT p
% so all constraints use current settings, and only borrow z_opt as z0.
%% =========================================================
warm_start_file = 'hzd_result_2026-07-20_12-06-36.mat';   % '' = cold start from x0
if ~isempty(warm_start_file)
    % Resolve relative to the repo's Results/ folder, independent of cwd
    % (same convention as the SAVE block below).
    ws_path = fullfile(fileparts(mfilename('fullpath')), '..', 'Results', warm_start_file);
    S  = load(ws_path);
    z0 = S.z_opt;
    fprintf('Warm start from %s\n', ws_path);
    fprintf('--- starting gait evaluated under CURRENT p (check its speed!) ---\n');
    inspect_solution(z0, p);
end

%% =========================================================
% OPTIMIZATION
%% =========================================================
[z_opt, fval, exitflag] = fmincon(...
    @(z) hzd_cost(z, p), ...
    z0, [], [], [], [], lb, ub, ...
    @(z) hzd_constraints(z, p), ...
    options);

report = inspect_solution(z_opt, p);

%% =========================================================
% SAVE RESULT
%% =========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
ts = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
result_file = fullfile(results_dir, ['hzd_result_', ts, '.mat']);
save(result_file, 'z_opt', 'p', 'fval', 'exitflag', 'report');
fprintf('Result saved to: %s\n', result_file);

animate_hzd_result(z_opt, p, 10);
