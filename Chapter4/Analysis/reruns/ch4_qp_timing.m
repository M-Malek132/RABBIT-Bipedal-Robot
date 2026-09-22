function T = ch4_qp_timing()
%CH4_QP_TIMING  Wall-clock cost of one control evaluation, per law.
%
%   T = ch4_qp_timing()
%
% WHY THIS EXISTS. Chapter 4's central thesis is that sampled-data effects
% dominate: the L1 estimator outruns the 1 kHz loop, and the cleanest fix
% measured anywhere in the chapter is a shorter control period. That
% recommendation is only meaningful if the loop can be closed at that rate,
% so the chapter has to report what one evaluation costs. It did not.
%
% WHAT IS MEASURED. One ch4_control call -- the QP solve plus the I/O
% linearization around it -- at states taken from ALONG one step of the
% reference gait, so the timing covers the problems the loop actually sees
% rather than one lucky state. Reported as median and 99th percentile: a 1 ms
% budget is missed by the tail, not by the median.
%
% CAVEAT, AND IT IS A REAL ONE. This is interpreted MATLAB calling quadprog on
% a laptop. It is an upper bound on a compiled implementation by a wide and
% unmeasured margin, and it is NOT a claim about embedded feasibility. What it
% does establish is (a) the RATIO between the laws, which is far less
% implementation dependent, and (b) whether the constrained laws cost
% materially more than the unconstrained one.

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
odir = fullfile(root, 'Results', 'reruns', 'ch4');
if ~exist(odir, 'dir'), mkdir(odir); end

[x0, alpha, p0] = ch4_load_gait();

% States along one step, from a baseline rollout at the nominal scale. Using a
% real rollout rather than the collocation nodes means the sampled states are
% the ones the closed loop visits, including the post-impact transient.
pb = ch4_run_params(p0, 'clfqp', 1.0, []);
st = ch4_step(x0, [], alpha, pb);
n_state = min(200, size(st.x, 2));
X  = st.x(:, round(linspace(1, size(st.x, 2), n_state)));
fprintf('timing over %d states along one step (%d repeats each)\n', size(X,2), 25);

laws  = {'clfqp', 'clfqp_con', 'rclfqp_con', 'l1', 'l1_con'};
n_rep = 25;
T = struct([]);

for k = 1:numel(laws)
    % mass_scale 1.5 so the constrained laws are actually pushing on their
    % rows; a law whose constraints are all slack solves a cheaper problem.
    p  = ch4_run_params(p0, laws{k}, 1.5, []);
    if ch4_is_stateful(p)
        eta0 = zeros(2 * p.ny, 1);
        xi = ch4_l1_state('init', p, eta0);
    else
        xi = [];
    end

    dt_us = nan(1, size(X, 2));
    for i = 1:size(X, 2)
        x = X(:, i);
        try
            ch4_control(x, xi, alpha, p);               % warm the path
            t0 = tic;
            for r = 1:n_rep, ch4_control(x, xi, alpha, p); end
            dt_us(i) = 1e6 * toc(t0) / n_rep;
        catch ME
            fprintf('   state %d skipped: %s\n', i, ME.identifier);
        end
    end
    d = dt_us(~isnan(dt_us));
    T(k).law    = laws{k};
    T(k).n      = numel(d);
    T(k).med_us = median(d);
    T(k).p90_us = prctile(d, 90);
    T(k).p99_us = prctile(d, 99);
    T(k).max_us = max(d);
    fprintf('%-12s n=%3d  median %7.1f  p90 %7.1f  p99 %7.1f  max %7.1f us\n', ...
            laws{k}, T(k).n, T(k).med_us, T(k).p90_us, T(k).p99_us, T(k).max_us);
end

fprintf('\nBudget: 1 kHz = 1000 us per period; 2 kHz = 500 us.\n');
for k = 1:numel(T)
    fprintf('%-12s  p99 = %6.1f%% of 1 ms, %6.1f%% of 0.5 ms   |   max = %6.1f%% of 1 ms\n', ...
            T(k).law, 100*T(k).p99_us/1000, 100*T(k).p99_us/500, 100*T(k).max_us/1000);
end

save(fullfile(odir, 'qp_timing.mat'), 'T');
fprintf('CH4_QP_TIMING_DONE\n');
end
