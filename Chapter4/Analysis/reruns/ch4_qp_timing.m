function T = ch4_qp_timing()
%CH4_QP_TIMING  Wall-clock cost of one control evaluation, per law.
%
%   T = ch4_qp_timing()
%
% WHY THIS EXISTS. Chapter 4's central thesis is that sampled-data effects
% dominate: the L1 estimator outruns the 1 kHz loop, and the cleanest fix
% measured anywhere in the chapter is a 0.5 ms period. That recommendation is
% only meaningful if the loop can actually be closed at that rate, so the
% chapter has to report what one evaluation costs. It did not.
%
% WHAT IS MEASURED. One ch4_control call -- the QP solve plus everything
% around it -- at states taken from along the reference gait, so the timing
% is of the problems the loop actually sees rather than of one lucky state.
% Reported as median and 99th percentile, because a 1 ms budget is missed by
% the tail, not by the median.
%
% CAVEAT, AND IT IS A REAL ONE. This is interpreted MATLAB calling quadprog
% on a laptop. It is an upper bound on a compiled implementation by a wide
% and unmeasured margin, and it is NOT a claim about embedded feasibility.
% What it does establish is the RATIO between the laws, which is
% implementation independent, and whether the margin at 0.5 ms is comfortable
% or nonexistent in this environment.

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
odir = fullfile(root, 'Results', 'reruns', 'ch4');
if ~exist(odir, 'dir'), mkdir(odir); end

g = ch4_load_gait();
laws = {'clfqp', 'clfqp_con', 'rclfqp_con', 'l1', 'l1_con'};
n_state = 40;      % states along the gait
n_rep   = 25;      % repeats per state

T = struct([]);
for k = 1:numel(laws)
    p = ch4_params('controller', laws{k}, 'mass_scale', 1.5);
    p = ch4_run_params(p);
    xi = ch4_l1_state(p);
    idx = round(linspace(1, size(g.X, 2), n_state));
    dt_us = zeros(1, numel(idx));
    for i = 1:numel(idx)
        x = g.X(:, idx(i));
        ch4_control(x, xi, g.alpha, p);                 % warm the path
        t0 = tic;
        for r = 1:n_rep, ch4_control(x, xi, g.alpha, p); end
        dt_us(i) = 1e6 * toc(t0) / n_rep;
    end
    T(k).law    = laws{k};
    T(k).med_us = median(dt_us);
    T(k).p99_us = prctile(dt_us, 99);
    T(k).max_us = max(dt_us);
    fprintf('%-12s median %7.1f us   p99 %7.1f us   max %7.1f us\n', ...
            laws{k}, T(k).med_us, T(k).p99_us, T(k).max_us);
end

fprintf('\nBudget: 1 kHz = 1000 us per period, 2 kHz = 500 us.\n');
for k = 1:numel(T)
    fprintf('%-12s  %5.1f%% of a 1 ms period,  %5.1f%% of a 0.5 ms period (p99)\n', ...
            T(k).law, 100*T(k).p99_us/1000, 100*T(k).p99_us/500);
end

save(fullfile(odir, 'qp_timing.mat'), 'T');
fprintf('CH4_QP_TIMING_DONE\n');
end
