function col_plot_grf(files, labels, outpng)
%COL_PLOT_GRF  Ground-reaction-force audit of one or more collocation gaits.
%
%   col_plot_grf(files, labels, outpng)
%
%   files  : char or cellstr of saved Results/col_result_*.mat gaits.
%   labels : cellstr of legend names, one per file.
%   outpng : output name written into Results/ (default 'col_grf.png').
%
% Example:
%   col_plot_grf({'col_result_A.mat','col_result_B.mat'}, ...
%                {'plain','friction+torque'}, 'col_grf_compare.png')
%
% Panels: Fz(s), Fx(s), the friction ratio |Fx|/Fz against mu_max on a log
% axis, and the friction cone in the (Fz,Fx) plane. The cone plot is the
% decisive one: any node outside the wedge means the stance foot slips, so
% the gait is not realizable on ground with that friction coefficient.
%
% Also prints a per-gait constraint audit (which Table-3.1 toggles were on,
% peak torque vs u_max, worst friction ratio vs mu_max, node counts).

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(root));

    if nargin < 3 || isempty(outpng), outpng = 'col_grf.png'; end
    if ischar(files), files = {files}; end
    if nargin < 2 || isempty(labels)
        labels = cellfun(@(f) f, files, 'UniformOutput', false);
    end

    n = numel(files);
    C = lines(max(n,3));
    f = figure('Color','w','Position',[50 50 1500 900]);

    mu = [];  Fzmax = 0;
    D = cell(1,n);
    for j = 1:n
        fn = files{j};
        if ~isfile(fn), fn = fullfile(root,'Results',fn); end
        S = load(fn);
        p = S.p;
        [X, U, Lam, ~, ~] = col_unpack(S.z_opt, p);
        N = p.N;  nq = p.nq;
        th = zeros(1,N);
        for k = 1:N, th(k) = theta_of_q(X(1:nq,k)); end
        s = (th - th(1)) / (p.theta_plus - th(1));

        d.s = s;  d.Lam = Lam;  d.U = U;
        d.ratio = abs(Lam(1,:)) ./ max(Lam(2,:), 1e-9);
        D{j} = d;
        mu = p.mu_max;
        Fzmax = max(Fzmax, max(Lam(2,:)));

        onf = @(b) char(string(b));
        fprintf('\n--- %s\n    %s\n', labels{j}, fn);
        fprintf('    toggles: torque=%s friction=%s grf=%s impulse=%s\n', ...
            onf(isfield(p,'enforce_torque')  && p.enforce_torque), ...
            onf(isfield(p,'enforce_friction')&& p.enforce_friction), ...
            onf(isfield(p,'enforce_grf')     && p.enforce_grf), ...
            onf(isfield(p,'enforce_impulse') && p.enforce_impulse));
        fprintf('    exitflag %d, cost %.4g\n', S.exitflag, S.fval);
        fprintf('    Fz  min/max      : %8.1f / %8.1f N\n', min(Lam(2,:)), max(Lam(2,:)));
        fprintf('    Fx  min/max      : %8.1f / %8.1f N\n', min(Lam(1,:)), max(Lam(1,:)));
        fprintf('    |Fx/Fz| max      : %8.3f   (mu_max = %.2f)  -> %s\n', ...
            max(d.ratio), mu, ternary(max(d.ratio) <= mu + 1e-6, 'OK', 'VIOLATED'));
        fprintf('    nodes outside cone: %d / %d\n', sum(d.ratio > mu + 1e-6), N);
        fprintf('    peak |torque|    : %8.1f Nm  (u_max = %.0f)  -> %s\n', ...
            max(abs(U(:))), p.u_max(1), ternary(max(abs(U(:))) <= p.u_max(1) + 1e-6, 'OK', 'VIOLATED'));
    end

    W = 32*9.81;   % RABBIT total weight [N], for scale on the Fz panel

    % --- 1) normal force -------------------------------------------------
    subplot(2,2,1); hold on; grid on;
    for j = 1:n
        plot(D{j}.s, D{j}.Lam(2,:), '-o','Color',C(j,:),'LineWidth',1.6,'MarkerSize',4);
    end
    yline(0,'k-','LineWidth',1.2);
    yline(W,'k--','body weight','LabelHorizontalAlignment','left');
    xlabel('phase s'); ylabel('F_z (N)');
    title('Normal ground reaction force');
    legend(labels,'Location','best'); legend boxoff;

    % --- 2) tangential force ---------------------------------------------
    subplot(2,2,2); hold on; grid on;
    for j = 1:n
        plot(D{j}.s, D{j}.Lam(1,:), '-o','Color',C(j,:),'LineWidth',1.6,'MarkerSize',4);
    end
    yline(0,'k-','LineWidth',1.2);
    xlabel('phase s'); ylabel('F_x (N)');
    title('Tangential (friction) force');

    % --- 3) friction ratio vs the cone limit -----------------------------
    subplot(2,2,3); hold on; grid on;
    for j = 1:n
        plot(D{j}.s, D{j}.ratio, '-o','Color',C(j,:),'LineWidth',1.6,'MarkerSize',4);
    end
    yline(mu,'r-','LineWidth',2);
    text(0.02, mu, sprintf('  \\mu_{max} = %.2f', mu), 'Color','r', ...
         'VerticalAlignment','bottom','FontWeight','bold');
    set(gca,'YScale','log');
    xlabel('phase s'); ylabel('|F_x| / F_z');
    title('Friction requirement (log scale) - above the red line the foot SLIPS');

    % --- 4) friction cone in the (Fz,Fx) plane ---------------------------
    subplot(2,2,4); hold on; grid on;
    Fzr = [0, Fzmax*1.05];
    patch([Fzr(1) Fzr(2) Fzr(2) Fzr(1)], [0 mu*Fzr(2) -mu*Fzr(2) 0], ...
          [0.85 0.95 0.85], 'EdgeColor','none');
    plot(Fzr,  mu*Fzr, 'r-','LineWidth',2);
    plot(Fzr, -mu*Fzr, 'r-','LineWidth',2);
    for j = 1:n
        in = D{j}.ratio <= mu + 1e-6;
        plot(D{j}.Lam(2, in), D{j}.Lam(1, in), 'o', 'Color',C(j,:), ...
             'MarkerFaceColor',C(j,:),'MarkerSize',6);
        plot(D{j}.Lam(2,~in), D{j}.Lam(1,~in), 'x', 'Color',C(j,:), ...
             'LineWidth',2,'MarkerSize',10);
    end
    xline(0,'k-'); yline(0,'k-');
    xlabel('F_z (N)'); ylabel('F_x (N)');
    title(sprintf('Friction cone \\mu=%.2f  (o inside = OK, x outside = slips)', mu));

    sgtitle('Direct collocation - stance-foot ground reaction forces', 'FontWeight','bold');
    exportgraphics(f, fullfile(root,'Results',outpng), 'Resolution', 130);
    fprintf('\nWrote Results/%s\n', outpng);
end

function o = ternary(cond, a, b)
    if cond, o = a; else, o = b; end
end
