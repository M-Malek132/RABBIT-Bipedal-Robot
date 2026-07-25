function R = ch3_report(z, p, opts)
%CH3_REPORT  Full diagnostic read-out for a collocation result.
%
%   R = ch3_report(z, p)
%   R = ch3_report(z, p, struct('stability', true, 'simulate', 3))
%
% This is the first thing to read after a solve. It answers, in order:
%   * did the transcription converge (residuals)
%   * what gait is it (speed, step length, duration, phase sweep)
%   * is it PHYSICALLY REALIZABLE (all four Table 3.1 quantities, MEASURED and
%     compared against the limits whether or not those limits were enforced)
%   * does it stay on the zero dynamics surface across the step
%   * is it actually STABLE (Poincare rho) -- optional, it costs 26 sims
%   * does it survive being simulated forward under the real controller
%
% MEASURE BEFORE YOU CONSTRAIN.  Every Table 3.1 quantity is reported with its
% measured value AND the limit, with a marker showing whether the limit is
% currently enforced. The thesis limits are ATRIAS numbers (63 kg, 50:1 gears,
% so its 5 Nm is 250 Nm at the joint); RABBIT is ~30 kg direct drive. Copying
% them across without measuring first produces an infeasible problem and a
% solver that fails for reasons that look like bugs.
%
% Inputs
%   z    : collocation decision vector
%   p    : parameter struct
%   opts : optional struct
%            .stability  (false) compute the Poincare spectral radius
%            .simulate   (0)     simulate this many steps under p.controller
%            .verbose    (true)  print
%
% Output
%   R : struct of everything printed, for programmatic use
%
% See also CH3_COL_EVAL, CH3_POINCARE, CH3_FORCES.

p = ch3_upgrade_params(p);   % saved results may predate fields read below

if nargin < 3, opts = struct(); end
if ~isfield(opts,'stability'), opts.stability = false; end
if ~isfield(opts,'simulate'),  opts.simulate  = 0;     end
if ~isfield(opts,'verbose'),   opts.verbose   = true;  end

E = ch3_col_eval(z, p);
[c, ceq] = ch3_col_constraints(z, p);
[X, T, alpha] = ch3_col_unpack(z, p);
N = E.N;

R = struct();
R.alpha = alpha; R.T = T; R.X = X; R.N = N;
R.max_ceq = max(abs(ceq));
R.max_c   = max(c);

% --- residuals by block ---------------------------------------------------
i = 0;
blk = {'node-1 (on Z, phase, contact)', 13; ...
       'Hermite-Simpson defects',       p.nx*(N-1); ...
       'node-N (on guard S)',           2; ...
       'periodicity through Delta',     13; ...
       'NEC1 walking rate',             1};
R.blocks = struct();
for b = 1:size(blk,1)
    n = blk{b,2};
    R.blocks.(matlab.lang.makeValidName(blk{b,1})) = max(abs(ceq(i+1:i+n)));
    i = i + n;
end

% --- gait -----------------------------------------------------------------
R.L_step = E.L_step;
R.speed  = E.L_step / E.T;
R.theta_sweep = [p.c_theta*X(1:p.nq,1), p.c_theta*X(1:p.nq,N)];
R.cost = ch3_col_cost(z, p);
R.int_u2 = E.int_u2;

% --- Table 3.1 quantities, measured --------------------------------------
lam_all = [E.lam, E.lamm];
R.torque_max  = max([abs(E.u(:)); abs(E.um(:))]);
R.Fz_min      = min(lam_all(2,:));
R.Fz_max      = max(lam_all(2,:));
pos = lam_all(2,:) > 1e-6;
if any(pos)
    R.mu_max = max(abs(lam_all(1,pos)) ./ lam_all(2,pos));
else
    R.mu_max = Inf;
end
R.impulse     = norm(E.impulse);
R.clearance   = E.sw_h(max(2, min(N-1, round((N+1)/2))));
R.sw_h_min    = min(E.sw_h(2:N-1));

% --- zero dynamics drift across the step ---------------------------------
eta_max = 0;
for k = 1:N
    [yk, ydk] = ch3_outputs(X(:,k), alpha, p);
    eta_max = max(eta_max, norm([yk; ydk], inf));
end
R.eta_max = eta_max;

% --- is this an actual trajectory, or a spurious discrete solution? -------
% Checked BEFORE anything else is believed: if the mesh is too coarse, every
% gait number above is fiction.
R.verify = ch3_col_verify(z, p, false);

% --- stance foot drift ----------------------------------------------------
foot0 = P_st(X(1:p.nq,1));
drift = 0;
for k = 1:N
    fk = P_st(X(1:p.nq,k));
    drift = max(drift, norm(fk - foot0, inf));
end
R.stance_drift = drift;

% --- stability ------------------------------------------------------------
R.rho = NaN;
if opts.stability
    [R.rho, ~, pinfo] = ch3_poincare(X(:,1), alpha, p);
    R.poincare = pinfo;
end

% --- forward simulation under the REAL controller -------------------------
R.sim = [];
if opts.simulate > 0
    R.sim = ch3_simulate(X(:,1), alpha, p, opts.simulate);
end

if ~opts.verbose, return; end

%% ------------------------------------------------------------- printing
line = repmat('-', 1, 74);
fprintf('\n%s\n CHAPTER 3 GAIT REPORT   (controller: %s, basis: %s, N = %d)\n%s\n', ...
        line, p.controller, p.basis, N, line);

fprintf(' CONVERGENCE\n');
fprintf('   max |ceq|                    %.3e\n', R.max_ceq);
fprintf('   max c (inequality)           %.3e\n', R.max_c);
fn = fieldnames(R.blocks);
for k = 1:numel(fn)
    fprintf('     %-30s %.3e\n', blk{k,1}, R.blocks.(fn{k}));
end

fprintf('\n GAIT\n');
fprintf('   step length                  %.4f m\n', R.L_step);
fprintf('   step duration                %.4f s\n', R.T);
fprintf('   walking speed                %.4f m/s   (v_des %.4f)\n', R.speed, p.v_des);
fprintf('   theta sweep                  %.4f -> %.4f   (target %.4f -> %.4f)\n', ...
        R.theta_sweep(1), R.theta_sweep(2), p.theta_minus, p.theta_plus);
fprintf('   cost  int||u||^2 / L_step    %.4f      (int||u||^2 = %.4f)\n', R.cost, R.int_u2);

fprintf('\n PHYSICAL REALIZABILITY  (Table 3.1)   [E] = enforced this solve\n');
prow('peak |torque|',      R.torque_max, p.limits.u_max,       '<=', 'Nm', p.limits.enable.torque);
prow('impact impulse',     R.impulse,    p.limits.impulse_max, '<=', 'Ns', p.limits.enable.impulse);
prow('friction demanded',  R.mu_max,     p.limits.mu_s,        '<=', '-',  p.limits.enable.friction);
prow('min vertical GRF',   R.Fz_min,     p.limits.Fz_min,      '>=', 'N',  p.limits.enable.grf);
prow('mid-step clearance', R.clearance,  p.limits.clearance,   '>=', 'm',  p.limits.enable.clearance);
fprintf('   (max Fz %.1f N, min swing-foot height over the step %.4f m)\n', R.Fz_max, R.sw_h_min);

fprintf('\n TRAJECTORY VALIDITY  (nodes vs a tight rollout -- read this FIRST)\n');
fprintf('   max |X_node - X_true|        %.3e   (tol %.1e)\n', R.verify.max_dev, p.verify_tol);
if R.verify.ok
    fprintf('   -> the nodes lie on a genuine trajectory; the numbers above are real.\n');
else
    fprintf('   -> *** REJECT: SPURIOUS DISCRETE SOLUTION ***\n');
    if ~R.verify.integrated
        fprintf('      The reference rollout could not even complete the step.\n');
    end
    fprintf('      The defects are satisfied but the nodes do not lie on any real\n');
    fprintf('      trajectory, so EVERY gait and Table 3.1 number above is fiction.\n');
    fprintf('      Cause is an under-resolved mesh. Re-solve with larger p.N_nodes\n');
    fprintf('      (currently %d).\n', N);
end

fprintf('\n INVARIANTS\n');
fprintf('   max |eta| over nodes         %.3e   (gait on the zero dynamics surface)\n', R.eta_max);
if R.eta_max > 1e-3
    fprintf('     ^ should be ~0. Large |eta| with a SMALL node-1 residual means the\n');
    fprintf('       mesh is too coarse, not that the invariance is wrong: u_ff makes\n');
    fprintf('       ydd = 0 exactly, so eta can only drift through discretization.\n');
end
fprintf('   stance-foot drift            %.3e m (pinned at every node)\n', R.stance_drift);

if opts.stability
    fprintf('\n ORBITAL STABILITY\n');
    fprintf('   Poincare rho                 %.4f   -> %s\n', R.rho, ...
            ternary(R.rho < 1, 'STABLE (attracting orbit)', 'UNSTABLE (periodic but falls)'));
    fprintf('   fixed-point residual         %.3e\n', R.poincare.fixed_point_residual);
    if R.rho >= 1
        fprintf('   NOTE: periodicity is satisfied but the orbit REPELS. Periodic is not\n');
        fprintf('         stable -- simulate it and the velocities grow step over step.\n');
    end
end

if opts.simulate > 0
    fprintf('\n FORWARD SIMULATION under "%s"  (%d steps requested)\n', p.controller, opts.simulate);
    fprintf('   steps completed              %d\n', R.sim.n_ok);
    fprintf('   outcome                      %s\n', R.sim.reason);
    if R.sim.n_ok > 0
        Ls = arrayfun(@(s) s.L_step, R.sim.steps);
        Ts = arrayfun(@(s) s.T,      R.sim.steps);
        fprintf('   step lengths                 %s\n', num2str(Ls, '%.3f '));
        fprintf('   step durations               %s\n', num2str(Ts, '%.3f '));
    end
end

fprintf('%s\n\n', line);

end

% ---------------------------------------------------------------------------
function prow(name, val, lim, sense, unit, enabled)
if strcmp(sense, '<=')
    ok = val <= lim;
else
    ok = val >= lim;
end
mark = ternary(ok, ' ok ', 'OVER');
if strcmp(sense,'>='), mark = ternary(ok, ' ok ', 'UNDER'); end
tag = ternary(enabled, '[E]', '[ ]');
fprintf('   %s %-22s %10.4f %s %-9.4f %-3s  %s\n', tag, name, val, sense, lim, unit, mark);
end

function o = ternary(cond, a, b)
if cond, o = a; else, o = b; end
end
