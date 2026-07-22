function R = col_crosscheck(col_file, controller, do_plot)
%COL_CROSSCHECK  Is a collocation gait realizable in closed loop?
%
%   R = col_crosscheck(col_file, controller, do_plot)
%
% A direct-collocation solution is an OPEN-LOOP trajectory: it satisfies the
% dynamics as defect constraints, but nothing guarantees the SHOOTING simulator
% under a feedback law reproduces it. This extracts the spline coeffs and the
% node-1 state from a saved col_result_*.mat, simulates ONE step under the
% chosen controller (p.controller = 'pd' or 'clfqp'), and compares the closed-
% loop trajectory to the collocation trajectory.
%
% Args:
%   col_file    Results/col_result_*.mat (default: newest). Must hold z_opt,p.
%   controller  'pd' (default) or 'clfqp'
%   do_plot     true to draw the overlay (default false; off under -batch)
%
% Returns R with the headline mismatches (SI units, rad):
%   R.q_rms, R.q_max        joint-position error, closed-loop vs collocation
%   R.Lstep_col/_sim, R.speed_col/_sim, R.Tstep_col/_sim
%   R.period_defect         ||reset(impact(x_end_sim)) - x_start(2:14)||inf
%   R.status                shooting-sim status (1 = valid step)

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(root));

    if nargin < 3 || isempty(do_plot),    do_plot = false;   end
    if nargin < 2 || isempty(controller), controller = 'pd'; end
    if nargin < 1 || isempty(col_file)
        d = dir(fullfile(root,'Results','col_result_*.mat'));
        assert(~isempty(d), 'No Results/col_result_*.mat; run rabbit_hzd_collocation first.');
        [~,i] = max([d.datenum]);
        col_file = fullfile(d(i).folder, d(i).name);
    end
    fprintf('Collocation gait : %s\ncontroller       : %s\n', col_file, controller);

    S = load(col_file);
    p = S.p;
    [X, ~, ~, coeffs, T] = col_unpack(S.z_opt, p);
    nq = p.nq;   N = p.N;

    % --- collocation reference trajectory (node states) ------------------
    tcol = linspace(0, T, N);
    Qcol = X(1:nq, :);                         % nq x N

    x_start = X(:,1);                          % node-1 state = closed-loop x0

    % --- closed-loop shooting sim under the chosen controller ------------
    % Own bounded integration (looser tol, capped horizon, try/catch) rather
    % than simulate_hzd_gait_full: the CLF-QP law solves a QP per ODE step, so
    % at the default 1e-6 tolerances a single step can crawl for tens of minutes
    % through the near-singular impact. Both controllers use these identical
    % settings so the comparison is fair.
    p.controller = controller;
    [tsim, xsim, status] = sim_step(coeffs, x_start, p);
    if status ~= 1
        warning('col_crosscheck:sim', ...
            '%s closed-loop step did NOT complete a valid impact (status=%d).', ...
            controller, status);
    end

    % --- resample sim to the collocation phase grid for comparison -------
    % Compare on a common NORMALIZED time grid tau in [0,1] (the two steps can
    % have slightly different durations). Guard against a degenerate sim.
    tau    = linspace(0,1,100);
    Qcol_i = interp1(tcol/T, Qcol.', tau).';                       % nq x 100
    if numel(tsim) >= 2
        Qsim_i  = interp1(tsim/tsim(end), xsim(:,1:nq), tau).';    % nq x 100
        err     = Qsim_i - Qcol_i;
        R.q_rms = sqrt(mean(err(:).^2));
        R.q_max = max(abs(err(:)));
    else
        Qsim_i  = nan(nq, numel(tau));
        R.q_rms = NaN;  R.q_max = NaN;
    end

    % --- gait metrics: collocation vs sim --------------------------------
    R.Tstep_col = T;   R.Tstep_sim = tsim(end);
    qN_col = X(1:nq,N);  swc = P_sw(qN_col); stc = P_st(qN_col);
    R.Lstep_col = swc(1) - stc(1);
    qN_sim = xsim(end,1:nq).';  sws = P_sw(qN_sim); sts = P_st(qN_sim);
    R.Lstep_sim = sws(1) - sts(1);
    R.speed_col = R.Lstep_col / R.Tstep_col;
    R.speed_sim = R.Lstep_sim / R.Tstep_sim;

    % --- how periodic is the closed-loop step? ---------------------------
    xend  = xsim(end,:).';
    xplus = rabbit_reset_map(rabbit_impact_map(xend));
    R.period_defect = max(abs(xplus(2:14) - x_start(2:14)));
    R.status = status;

    % --- report -----------------------------------------------------------
    fprintf('\nclosed-loop vs collocation:\n');
    fprintf('  status              : %d (1 = valid step)\n', R.status);
    fprintf('  joint pos err  rms  : %.4f rad\n', R.q_rms);
    fprintf('  joint pos err  max  : %.4f rad\n', R.q_max);
    fprintf('  step length col/sim : %.4f / %.4f m\n', R.Lstep_col, R.Lstep_sim);
    fprintf('  step time   col/sim : %.4f / %.4f s\n', R.Tstep_col, R.Tstep_sim);
    fprintf('  speed       col/sim : %.4f / %.4f m/s\n', R.speed_col, R.speed_sim);
    fprintf('  closed-loop periodicity defect (inf) : %.4e\n', R.period_defect);

    if do_plot
        figure('Name','col vs closed-loop');
        names = {'q1 stance-hip','q2 stance-knee','q3 swing-hip','q4 swing-knee'};
        for j = 1:4
            subplot(2,2,j); hold on; grid on;
            plot(tau, Qcol_i(3+j,:), 'k-',  'DisplayName','collocation');
            plot(tau, Qsim_i(3+j,:), 'r--', 'DisplayName',[controller ' sim']);
            title(names{j}); xlabel('phase'); ylabel('rad');
            if j==1, legend('Location','best'); end
        end
    end
end

% -------------------------------------------------------------------------
function [t, x, status] = sim_step(coeffs, x_start, p)
%SIM_STEP  One closed-loop step under p.controller, with BOUNDED cost.
% Looser tolerances + a MaxStep ceiling + a capped horizon keep the CLF-QP
% (QP-per-step) sim tractable, and a try/catch means a divergent closed loop
% returns status=-1 instead of hanging. Same impact event as the main sims.
    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;
    Tcap = min(p.T_max, 1.5);

    opts = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',0.02, ...
                  'Events',@impact_event_wrapper);
    xi0 = [x_start; 0];
    try
        sol = ode45(@(t,xi) hzd_ode_rhs(t,xi,coeffs,theta_minus,theta_plus,p), ...
                    [0 Tcap], xi0, opts);
    catch ME
        warning('col_crosscheck:odefail','ode45 failed: %s', ME.message);
        t = 0; x = x_start(:).'; status = -1; return;
    end

    t = sol.x(:);
    x = sol.y(1:2*nq, :).';
    impacted = isfield(sol,'ie') && ~isempty(sol.ie);
    status = (impacted && ~any(isnan(x(end,:)))) * 1 + ...
             (~(impacted && ~any(isnan(x(end,:))))) * -1;
end
