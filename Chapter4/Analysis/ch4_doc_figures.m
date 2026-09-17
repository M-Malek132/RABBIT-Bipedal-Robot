function ch4_doc_figures()
%CH4_DOC_FIGURES  Redraw the Chapter-4 report figures into docs/figures.
%
%   ch4_doc_figures()
%
% The report (docs/ch4_report_fa.tex) shows eight figures drawn from four saved
% result sets. This redraws them from the .mat files, so a change of figure
% style does not mean re-running the sweeps -- which take the better part of an
% hour, and which MATLAB here rarely survives in one session. Then rebuild:
%
%   ch4_doc_figures
%   !cd docs && latexmk -xelatex ch4_report_fa.tex
%
%   docs/figures/          result set (Results/)               drawn by, fig
%   ch4_robust_clf.png     ch4_result_2026-09-13_20-45-50      ch4_plot_uncertainty 1
%   ch4_robust_torso.png     "                                  ch4_plot_uncertainty 4
%   ch4_case4_clf.png      ch4_result_2026-09-14_13-50-53      ch4_plot_uncertainty 1
%   ch4_l1_clf.png         ch4_result_2026-09-14_22-33-20      ch4_plot_uncertainty 1
%   ch4_l1_theta.png         "                                  ch4_plot_uncertainty 5
%   ch4_load_random.png    ch4_result_2026-09-14_22-37-57      ch4_plot_load 1
%   ch4_load_torque.png      "                                  ch4_plot_load 2
%   ch4_load_torso.png       "                                  ch4_plot_load 3
%
% THE RESULT SETS ARE PINNED, NOT THE NEWEST. The report's numbers were read
% from these runs, and its figures must show the same runs. After a rerun,
% update the stamps here together with the numbers in the report.
%
% The widths the figures are drawn at (ch4_figure) are the widths the .tex sets
% them at: a column for the CLF, estimator and robust phase-portrait figures,
% the page for the three load figures.
%
% See also CH4_PLOT_UNCERTAINTY, CH4_PLOT_LOAD, CH4_FIGURE, CH4_MAIN.

root   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
res    = @(stamp) fullfile(root, 'Results', sprintf('ch4_result_%s.mat', stamp));
figdir = fullfile(root, 'docs', 'figures');

R = load(res('2026-09-13_20-45-50'), 'robust', 'p');
figs = ch4_plot_uncertainty(R.robust, R.p);
save_as(figs(1), 'ch4_robust_clf.png');
save_as(figs(4), 'ch4_robust_torso.png');
close(figs);

R = load(res('2026-09-14_13-50-53'), 'case4', 'p_case4');
figs = ch4_plot_uncertainty(R.case4, R.p_case4, '', {'IV'});
save_as(figs(1), 'ch4_case4_clf.png');
close(figs);

R = load(res('2026-09-14_22-33-20'), 'l1', 'p');
figs = ch4_plot_uncertainty(R.l1, R.p);
save_as(figs(1), 'ch4_l1_clf.png');
save_as(figs(5), 'ch4_l1_theta.png');
close(figs);

% The nominal orbit under the load portraits is not saved with the study.
% ch4_main takes it from a two-step nominal rollout of the baseline; redo that.
R  = load(res('2026-09-14_22-37-57'), 'load', 'p', 'x0', 'alpha');
pb = R.p;
pb.controller        = 'clfqp';
pb.uncertainty       = struct('mass_scale', 1, 'load_mass', 0);
pb.load_random_range = [];
sim_b = ch4_simulate(R.x0, R.alpha, pb, 2);
figs = ch4_plot_load(R.load, R.p, '', struct('X_orbit', sim_b.x(:, 1:2:end)));
save_as(figs(1), 'ch4_load_random.png');
save_as(figs(2), 'ch4_load_torque.png');
save_as(figs(3), 'ch4_load_torso.png');
close(figs);

    function save_as(f, name)
        exportgraphics(f, fullfile(figdir, name), 'Resolution', 600);
        fprintf(' ch4_doc_figures: %s\n', name);
    end
end
