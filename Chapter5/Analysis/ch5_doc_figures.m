function ch5_doc_figures()
%CH5_DOC_FIGURES  Redraw the Chapter-5 report figures into docs/figures.
%
%   ch5_doc_figures()
%
% The report (docs/ch5_report_fa.tex) shows four figures from one result set.
% This redraws them from its .mat, so a change of figure style does not mean
% re-running the studies. Then rebuild:
%
%   ch5_doc_figures
%   !cd docs && latexmk -xelatex ch5_report_fa.tex
%
%   docs/figures/                 drawn by
%   ch5_springmass_fig53.png      ch5_plot_springmass, fig 1
%   ch5_springmass_barrier.png    ch5_plot_springmass, fig 2
%   ch5_pendulum_fig54.png        ch5_plot_pendulum, fig 1
%   ch5_pendulum_barrier.png      ch5_plot_pendulum, fig 2
%
% THE RESULT SET IS PINNED, NOT THE NEWEST: ch5_result_2026-08-17_14-04-43 is
% the run the report's numbers were read from, and its figures must show the
% same run. After a rerun, update the stamp here together with the report.
%
% The plot functions save under the names the report uses, so they write
% straight into docs/figures.
%
% See also CH5_PLOT_SPRINGMASS, CH5_PLOT_PENDULUM, CH5_FIGURE, CH5_MAIN.

root   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(root, 'docs', 'figures');

R = load(fullfile(root, 'Results', 'ch5_result_2026-08-17_14-04-43.mat'), ...
         'springmass', 'pendulum');

close(ch5_plot_springmass(R.springmass, figdir));
fprintf(' ch5_doc_figures: spring-mass figures written\n');
close(ch5_plot_pendulum(R.pendulum, figdir));
fprintf(' ch5_doc_figures: pendulum figures written\n');
end
