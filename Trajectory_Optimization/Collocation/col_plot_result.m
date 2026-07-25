function col_plot_result(col_file, outpng)
%COL_PLOT_RESULT  Visual + numeric report of a saved direct-collocation gait.
%
%   col_plot_result(col_file, outpng)
%
%   col_file : saved Results/col_result_*.mat (default: newest).
%   outpng   : output name written into Results/ (default 'col_report.png').
%
% Prints the feasibility/task summary (exitflag, max|ceq|, max c, speed, step
% length, peak torque, GRF range, friction ratio) and writes a six-panel
% figure: virtual constraints vs the B-spline, node torques against the
% torque box, stance contact force, the zero-dynamics phase portrait,
% swing-foot clearance, and a stick figure of the step.
%
% See also COL_PLOT_GRF for a dedicated ground-reaction-force / friction-cone
% audit across several gaits.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(root));

    if nargin < 1 || isempty(col_file)
        d = dir(fullfile(root,'Results','col_result_*.mat'));
        [~,i] = max([d.datenum]);
        col_file = fullfile(d(i).folder, d(i).name);
    elseif ~isfile(col_file)
        col_file = fullfile(root,'Results', col_file);
    end
    if nargin < 2 || isempty(outpng), outpng = 'col_report.png'; end

    S = load(col_file);
    p = S.p;
    [X, U, Lam, coeffs, T] = col_unpack(S.z_opt, p);
    nq = p.nq;  nu = p.nu;  N = p.N;

    % ---- feasibility ----------------------------------------------------
    [c, ceq] = col_constraints(S.z_opt, p);

    Q  = X(1:nq,:);   dQ = X(nq+1:end,:);
    th = zeros(1,N);  dth = zeros(1,N);
    for k = 1:N
        th(k)  = theta_of_q(Q(:,k));
        dth(k) = [0 0 1 1 0.5 0 0] * dQ(:,k);
    end
    s = (th - th(1)) / (p.theta_plus - th(1));

    swN = P_sw(Q(:,N));  stN = P_st(Q(:,N));
    Lstep = swN(1) - stN(1);

    % swing foot height over the step
    zsw = zeros(1,N);
    for k = 1:N, ps = P_sw(Q(:,k)); zsw(k) = ps(2); end

    fprintf('\n================ DIRECT COLLOCATION RESULT ================\n');
    fprintf('file            : %s\n', col_file);
    fprintf('nodes N         : %d   (decision vars = %d)\n', N, numel(S.z_opt));
    fprintf('exitflag        : %d\n', S.exitflag);
    fprintf('cost (int u^2/L): %.6g\n', S.fval);
    fprintf('max |ceq|       : %.3e   (equality feasibility)\n', max(abs(ceq)));
    fprintf('max c           : %.3e   (inequality, <=0 is satisfied)\n', max(c));
    fprintf('-----------------------------------------------------------\n');
    fprintf('step duration T : %.4f s\n', T);
    fprintf('step length     : %.4f m\n', Lstep);
    fprintf('avg speed       : %.4f m/s   (target %.4f)\n', Lstep/T, p.v_des);
    fprintf('theta  start/end: %.4f -> %.4f rad\n', th(1), th(N));
    fprintf('peak |torque|   : %s Nm\n', mat2str(round(max(abs(U),[],2)',1)));
    fprintf('peak swing clr  : %.4f m\n', max(zsw));
    fprintf('GRF Fz min/max  : %.1f / %.1f N\n', min(Lam(2,:)), max(Lam(2,:)));
    fprintf('max |Fx/Fz|     : %.3f   (friction ratio)\n', max(abs(Lam(1,:)./max(Lam(2,:),1e-6))));
    fprintf('===========================================================\n\n');

    % ---- figure ---------------------------------------------------------
    f = figure('Color','w','Position',[50 50 1500 950]);
    tl = {'stance hip q4','stance knee q5','swing hip q6','swing knee q7'};

    % 1) actuated joints: node values vs the B-spline they must lie on
    subplot(2,3,1); hold on; grid on;
    ss = linspace(0,1,200);
    C = lines(4);
    for i = 1:nu
        yd = arrayfun(@(x) bspline_eval(coeffs(i,:), x, p), ss);
        plot(ss, yd, '-', 'Color', C(i,:), 'LineWidth', 1.5);
        plot(s, Q(3+i,:), 'o', 'Color', C(i,:), 'MarkerSize', 4, ...
             'MarkerFaceColor', C(i,:));
    end
    xlabel('phase s'); ylabel('rad');
    title('Virtual constraints: nodes (o) on B-spline yd(s)');
    legend(tl, 'Location','best'); legend boxoff;

    % 2) torques
    subplot(2,3,2); hold on; grid on;
    for i = 1:nu, plot(s, U(i,:), '-o', 'Color', C(i,:), 'MarkerSize',3); end
    if isfield(p,'u_max') && ~isempty(p.u_max)
        yline( p.u_max(1),'k--'); yline(-p.u_max(1),'k--');
    end
    xlabel('phase s'); ylabel('N\cdotm'); title('Node torques u_k');

    % 3) contact forces
    subplot(2,3,3); hold on; grid on;
    plot(s, Lam(2,:), '-o','LineWidth',1.5,'MarkerSize',3);
    plot(s, Lam(1,:), '-s','LineWidth',1.5,'MarkerSize',3);
    yline(0,'k:');
    xlabel('phase s'); ylabel('N'); title('Stance contact force \lambda');
    legend({'F_z (normal)','F_x (tangential)'},'Location','best'); legend boxoff;

    % 4) phase portrait of the zero dynamics
    subplot(2,3,4); hold on; grid on;
    plot(th, dth, '-o','LineWidth',1.5,'MarkerSize',4);
    plot(th(1), dth(1), 'go','MarkerFaceColor','g','MarkerSize',9);
    plot(th(N), dth(N), 'rs','MarkerFaceColor','r','MarkerSize',9);
    xlabel('\theta (rad)'); ylabel('d\theta/dt (rad/s)');
    title('Zero-dynamics phase portrait (o start, \square impact)');

    % 5) swing foot height
    subplot(2,3,5); hold on; grid on;
    plot(s, zsw, '-o','LineWidth',1.5,'MarkerSize',3);
    yline(0,'k-','LineWidth',1.5);
    xlabel('phase s'); ylabel('z_{swing} (m)'); title('Swing-foot clearance');

    % 6) stick figure over the step
    subplot(2,3,6); hold on; grid on; axis equal;
    plot([-0.2 Lstep+0.4],[0 0],'k-','LineWidth',2);
    idx = round(linspace(1,N,7));
    g = linspace(0.82,0.0,numel(idx));
    for j = 1:numel(idx)
        q = Q(:,idx(j));
        b = get_body_points(q);
        col = [g(j) g(j) g(j)];
        plot([b.hip(1) b.stance_knee(1) b.stance_foot(1)], ...
             [b.hip(2) b.stance_knee(2) b.stance_foot(2)], '-', 'Color',col,'LineWidth',2);
        plot([b.hip(1) b.swing_knee(1)  b.swing_foot(1)], ...
             [b.hip(2) b.swing_knee(2)  b.swing_foot(2)], '-', 'Color',col,'LineWidth',2);
        plot([b.hip(1) b.torso_top(1)], [b.hip(2) b.torso_top(2)], '-','Color',col,'LineWidth',3);
    end
    xlabel('x (m)'); ylabel('z (m)'); title('Gait over one step (light\rightarrowdark)');

    sgtitle(sprintf('Direct collocation (Hermite-Simpson), N=%d  |  v=%.3f m/s, T=%.3f s, L=%.3f m  |  max|ceq|=%.1e', ...
        N, Lstep/T, T, Lstep, max(abs(ceq))), 'FontWeight','bold');

    exportgraphics(f, fullfile(root,'Results',outpng), 'Resolution', 130);
    fprintf('Wrote Results/%s\n', outpng);
end
