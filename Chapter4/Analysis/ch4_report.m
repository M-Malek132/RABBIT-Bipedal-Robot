function R = ch4_report(x0, alpha, p, opts)
%CH4_REPORT  Summarize one controller running against one perturbed model.
%
%   R = ch4_report(x0, alpha, p)
%   R = ch4_report(x0, alpha, p, opts)
%
% ch3_report answers "is this gait any good?". This answers a different
% question: "what happened when the controller's model was wrong?" -- so
% alongside the usual tracking, torque and contact numbers it reports the
% uncertainty the controller was actually facing, and, for the adaptive laws,
% what it made of it.
%
% THE THREE BLOCKS AND WHAT EACH IS FOR
%
%   TRACKING     max and final ||eta||, and the per-step decay ratio. The
%                chapter's claim is about the RATE being preserved under
%                perturbation, so the decay ratio -- not the raw error -- is
%                the number that tests it.
%
%   ACTUATION    peak torque, and how far past its own box the law went. A
%                nonzero excess is expected for 'l1_con' (saturation binds on
%                mu1 only, Section 4.2.3) and is a bug for 'rclfqp_con'.
%
%   CONTACT      Fz_min and the required friction, computed from the TRUE
%                model, plus grf_pred_error -- how far the controller's own
%                nominal prediction of the contact force was from the real one.
%                That last number is the concrete size of the gap Remark 4.4
%                warns about, and it is why constraints written on the nominal
%                model are not robust constraints.
%
%   UNCERTAINTY  the realized ||Delta1||, ||Delta2|| along the run, checked
%                against the bounds the robust controller was designed with.
%                If bound_violated is true the robust guarantee did not apply
%                to this run and its results should not be cited as evidence
%                for it -- the controller was outside its own hypothesis.
%
% Inputs
%   x0, alpha : the gait (ch4_load_gait)
%   p         : parameter struct, including p.controller and p.uncertainty
%   opts      : .n_steps (default 3) .verbose (default true) .stability
%               (default false; runs a Poincare check, which costs more steps)
%
% Output
%   R : struct of every number printed, plus .sim and .F for further analysis
%
% See also CH4_COMPARE_CONTROLLERS, CH4_FORCES, CH3_REPORT.

if nargin < 4, opts = struct(); end
opts = fill_defaults(opts, struct('n_steps', 3, 'verbose', true, ...
                                  'stability', false));

sim = ch4_simulate(x0, alpha, p, opts.n_steps);

R = struct();
R.controller = p.controller;
R.uncertainty = p.uncertainty;
R.n_steps_ok = sim.n_ok;
R.failed     = sim.failed;
R.reason     = sim.reason;
R.sim        = sim;

if sim.n_ok == 0
    if opts.verbose
        fprintf('\n  ch4_report: no completed steps -- %s\n\n', sim.reason);
    end
    R.F = [];
    return;
end

F = ch4_forces(sim.t, sim.x, alpha, p, sim.xi, sim.t_xi);
R.F = F;

eta_n = vecnorm([F.y; F.ydot], 2, 1);

R.max_eta   = max(eta_n);
R.final_eta = eta_n(end);

% per-step decay: ||eta|| at the end of each step over its value just after the
% preceding impact. The chapter's claim is that this stays put as the model is
% perturbed.
R.step_decay = zeros(1, sim.n_ok);
t_edge = 0;
for k = 1:sim.n_ok
    t_end = t_edge + sim.steps(k).T;
    in_k  = F.t >= t_edge - 1e-12 & F.t <= t_end + 1e-12;
    e_k   = eta_n(in_k);
    if numel(e_k) > 1
        R.step_decay(k) = e_k(end) / max(e_k(1), realmin);
    end
    t_edge = t_end;
end

R.peak_torque   = F.torque_max;
R.int_u2        = F.int_u2;
R.box           = box_for(p);
R.box_excess    = max(F.torque_max - R.box, 0);
R.delta_max     = F.delta_max;
R.qp_infeasible = F.qp_infeasible;

R.Fz_min         = F.Fz_min;
R.mu_required    = F.mu_max;
R.grf_pred_error = F.grf_pred_error;

R.Delta1_max = max(vecnorm(F.Delta1, 2, 1));
R.Delta2_max = max(F.Delta2_n);
R.bound_violated = (R.Delta1_max > p.rclf.delta1_max) || ...
                   (R.Delta2_max > p.rclf.delta2_max);

R.step_T      = [sim.steps.T];
R.step_L      = [sim.steps.L_step];
R.v_avg       = sum(R.step_L) / sum(R.step_T);
R.impulse_max = max(cellfun(@(a) norm(a), {sim.steps.impulse}));

if ch4_is_stateful(p)
    R.theta_hat_max = max(vecnorm(F.theta_hat, 2, 1));
    R.theta_err_max = F.theta_err_max;
    R.mu2_max       = max(vecnorm(F.mu2, 2, 1));
else
    R.theta_hat_max = NaN; R.theta_err_max = NaN; R.mu2_max = NaN;
end

if opts.stability
    R.stability = ch3_stability(sim.x_final, alpha, p);
end

if opts.verbose, print_report(R, p); end

end

% ---------------------------------------------------------------------------
function print_report(R, p)
w = 74;
fprintf('\n%s\n', repmat('=',1,w));
fprintf(' CH4 REPORT  controller "%s"  |  mass scale %.2f, load %.1f kg\n', ...
        R.controller, R.uncertainty.mass_scale, R.uncertainty.load_mass);
fprintf('%s\n', repmat('-',1,w));

fprintf(' steps completed        %d%s\n', R.n_steps_ok, ...
        tern(R.failed, sprintf('   (stopped: %s)', R.reason), ''));
fprintf(' mean step / speed      T = %.4f s,  L = %.4f m,  v = %.4f m/s\n', ...
        mean(R.step_T), mean(R.step_L), R.v_avg);

fprintf('\n TRACKING\n');
fprintf('   max ||eta||          %.4e\n', R.max_eta);
fprintf('   final ||eta||        %.4e\n', R.final_eta);
fprintf('   per-step decay       %s\n', mat2str(round(R.step_decay,4)));

fprintf('\n ACTUATION\n');
fprintf('   peak |u|             %.2f Nm   (box %.0f Nm)\n', R.peak_torque, R.box);
fprintf('   over box             %.2f Nm%s\n', R.box_excess, ...
        tern(R.box_excess > 1e-6 && any(strcmpi(p.controller,{'l1_con'})), ...
             '   <- expected: saturation binds on mu1 only (4.2.3)', ''));
fprintf('   int ||u||^2          %.4e\n', R.int_u2);
fprintf('   CLF slack max        %.4e\n', R.delta_max);
fprintf('   QP infeasible        %d samples\n', R.qp_infeasible);

fprintf('\n CONTACT  (TRUE model, not the controller''s prediction)\n');
fprintf('   min Fz               %.2f N%s\n', R.Fz_min, ...
        tern(R.Fz_min <= 0, '   <- foot left the ground', ''));
fprintf('   required friction    %.4f\n', R.mu_required);
fprintf('   max impact impulse   %.3f Ns\n', R.impulse_max);
fprintf('   nominal GRF error    %.2f N   <- size of the Remark 4.4 gap\n', ...
        R.grf_pred_error);

fprintf('\n UNCERTAINTY FACED\n');
fprintf('   max ||Delta1||       %8.2f   (bound %.1f)%s\n', ...
        R.Delta1_max, p.rclf.delta1_max, ...
        tern(R.Delta1_max > p.rclf.delta1_max, '  EXCEEDED', ''));
fprintf('   max ||Delta2||       %8.4f   (bound %.3f)%s\n', ...
        R.Delta2_max, p.rclf.delta2_max, ...
        tern(R.Delta2_max > p.rclf.delta2_max, '  EXCEEDED', ''));
if R.bound_violated && any(strncmpi(p.controller, {'rclfqp'}, 6))
    fprintf('   *** the run left the uncertainty set the controller assumed,\n');
    fprintf('   *** so the robust guarantee does not cover it.\n');
end

if ~isnan(R.theta_hat_max)
    fprintf('\n L1 ESTIMATOR\n');
    fprintf('   max ||theta_hat||    %.4f\n', R.theta_hat_max);
    fprintf('   max ||theta err||    %.4f   (vs the true Delta1 + Delta2 mu)\n', ...
            R.theta_err_max);
    fprintf('   max ||mu2||          %.4f   (the adaptive contribution)\n', ...
            R.mu2_max);
end

fprintf('%s\n\n', repmat('=',1,w));
end

% ---------------------------------------------------------------------------
function b = box_for(p)
if any(strcmpi(p.controller, {'l1', 'l1_con'}))
    b = p.l1.u_max;
else
    b = p.limits.u_max;
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function s = fill_defaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i}))
        s.(f{i}) = d.(f{i});
    end
end
end
