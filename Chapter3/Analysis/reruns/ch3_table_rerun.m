function ch3_table_rerun()
%CH3_TABLE_RERUN  The Chapter-3 report's controller table (tab:ctrl) on posture_195,
% rerun with physical validity. Same setup as the original run (2026-09-16):
% 1 kHz ZOH, the gait's own CARE CLF at eps 0.5, only clfqp_con told about the
% box, p.T_max = 0.5 (a step not at the guard within 1.0 s fails), 4 steps, and
% a step counts as a fall when it is shorter than p.step_len_min, its tracking
% error exceeds 1e3, or the state is not finite. Boxes: 60% of the gait's peak
% for all three laws, then 80% and 100% for clfqp_con. One checkpoint per
% (configuration, step), so a crashed session resumes.
%
% Output: Results/reruns/ch3/ (one .mat per step, table.log).
ROOT = [fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))) filesep];
SPD  = [fullfile(ROOT, 'Results', 'reruns', 'ch3') filesep];
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = [SPD 'table.log'];

S = load(fullfile(ROOT, 'Results', 'ch3_gait_posture_195.mat'));
z = S.z; p = ch3_upgrade_params(S.p);
E = ch3_col_eval(z, p);
[X, ~, alpha] = ch3_col_unpack(z, p);
peak_ff = max([abs(E.u(:)); abs(E.um(:))]);
p.control_dt = 1e-3;
p.T_max = 0.5;
n_steps = 4;

configs = {'iolin_pd', 0.6; 'clfqp', 0.6; 'clfqp_con', 0.6; 'clfqp_con', 0.8; 'clfqp_con', 1.0};
logline(logf, 'start %s | eps %.2f %s | peak ff %.1f Nm | mu_s %.2f', datestr(now), ...
        p.eps, p.clf_construction, peak_ff, p.limits.mu_s);
for c = 1:size(configs, 1)
    [name, frac] = configs{c, :};
    u_box = frac * peak_ff;
    pc = p; pc.controller = name; pc.limits.u_max = u_box;
    pc.limits.enable.torque = strcmp(name, 'clfqp_con');
    x = X(:, 1);
    valid_run = 0; still_valid = true;
    for k = 1:n_steps
        f = sprintf('%s%s_box%03.0f_%d.mat', SPD, name, 100*frac, k);
        if exist(f, 'file')
            L = load(f); r = L.r;
        else
            if ~all(isfinite(x)), break; end
            s = ch3_step(x, alpha, pc);
            r = struct('name', name, 'box', u_box, 'k', k, 'ok', s.ok, 'T', s.T, ...
                       'L', s.L_step, 'x_next', s.x_next, 'peak', NaN, 'eta', NaN, ...
                       'delta', NaN, 'qpfail', NaN);
            try
                F = ch3_forces(s.t, s.x, alpha, pc);
                r.peak = F.torque_max; r.eta = max(vecnorm([F.y; F.ydot], 2, 1));
                r.delta = F.delta_max; r.qpfail = F.qp_infeasible;
            catch err
                r.forces_error = err.message;
            end
            r.V = ch3_validity(struct('steps', s), pc);
            r.held_peak = max(abs(s.u(:)));
            r.fell = ~s.ok || s.L_step < p.step_len_min || ~all(isfinite(s.x_next)) || ~(r.eta < 1e3);
            save(f, 'r');
        end
        step_valid = r.V.valid_steps == 1;
        if still_valid && step_valid && ~r.fell, valid_run = valid_run + 1; else, still_valid = false; end
        logline(logf, ['%-9s box %5.1f step %d ok %d T %.4f L %.4f | peak %.1f | max|eta| %.3e | ' ...
                       'delta %.2e | QPfail %d | fell %d | step valid %d (%s) minFz %.1f maxmu %.3f bad %.2f%%'], ...
                name, u_box, k, r.ok, r.T, r.L, r.peak, r.eta, r.delta, r.qpfail, r.fell, ...
                step_valid, r.V.first_kind, r.V.Fz_min, r.V.mu_max, 100 * r.V.frac_invalid);
        x = r.x_next;
        if r.fell, break; end
    end
    logline(logf, 'SUMMARY %-9s box %5.1f: steps before fall %d, valid steps %d', name, u_box, ...
            k - double(r.fell), valid_run);
end
logline(logf, 'DONE');
fprintf('TABLE_DONE\n');
end

function logline(f, fmt, varargin)
fid = fopen(f, 'a'); fprintf(fid, [fmt '\n'], varargin{:}); fclose(fid);
end
