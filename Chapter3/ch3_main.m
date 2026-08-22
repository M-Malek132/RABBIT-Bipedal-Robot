function out = ch3_main(varargin)
%CH3_MAIN  Run the Chapter-3 pipeline end to end.
%
%   out = ch3_main()
%   out = ch3_main('v_des', 0.30, 'controller', 'clfqp_con', ...)
%
% Walks the eight stages in dependency order:
%
%   (1) hybrid model            ch3_control_affine / ch3_guard / ch3_impact
%   (2) virtual constraints     ch3_phase / ch3_bezier / ch3_outputs
%   (4) I/O linearization       ch3_io_lin              <- needed by (3)
%   (3) offline optimization    ch3_col_solve           -> alpha
%   (5) PD baseline             ch3_ctrl_pd
%   (6) RES-CLF                 ch3_res_clf
%   (7) CLF-QP                  ch3_ctrl_clf_qp
%   (8) constrained CLF-QP      ch3_ctrl_clf_qp(..., true)
%
% Note the ordering: stage 3 is solved AFTER stage 4 is available, because the
% collocation transcription uses u_ff -- the feedforward from the I/O
% linearization -- as its input. The chapter presents 3 before 4 for narrative
% reasons; the code cannot.
%
% Any ch3_params field can be overridden by name/value.
%
% Output
%   out : struct .p .z_opt .alpha .solve .report
%
% See also CH3_PARAMS, CH3_COL_SOLVE, CH3_REPORT, CH3_COMPARE_CONTROLLERS.

p = ch3_params(varargin{:});

results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
stamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS'); %#ok<TNOW1,DATST>

if isempty(p.checkpoint_file)
    p.checkpoint_file = fullfile(results_dir, sprintf('ch3_ckpt_%s.mat', stamp));
end

fprintf('\n================ CHAPTER 3 PIPELINE ================\n');
fprintf(' basis %s (degree %d), controller %s, N = %d nodes\n', ...
        p.basis, p.bez_deg, p.controller, p.N_nodes);
fprintf(' v_des %.3f m/s, eps %.2f, CLF via %s\n', p.v_des, p.eps, p.clf_construction);

%% --- stage 3: solve for alpha -------------------------------------------
% The solve runs under pure feedforward regardless of p.controller: the gait
% is designed ON the zero dynamics surface, where the feedback term is
% multiplying zero. p.controller selects what RUNS the resulting gait.
[z_opt, solve_out] = ch3_col_solve(p);

%% --- report --------------------------------------------------------------
R = ch3_report(z_opt, p, struct('stability', true, 'simulate', 5));

%% --- save ---------------------------------------------------------------
alpha = solve_out.alpha; %#ok<NASGU>
fname = fullfile(results_dir, sprintf('ch3_result_%s.mat', stamp));
save(fname, 'z_opt', 'solve_out', 'p', 'R');
fprintf(' Saved: %s\n', fname);

if exist(p.checkpoint_file, 'file'), delete(p.checkpoint_file); end

out = struct('p', p, 'z_opt', z_opt, 'alpha', solve_out.alpha, ...
             'solve', solve_out, 'report', R, 'file', fname);

end
