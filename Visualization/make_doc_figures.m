function make_doc_figures()
%MAKE_DOC_FIGURES  Regenerate the result figures embedded in docs/.
%
%   make_doc_figures()
%
% Writes four PNGs into docs/figures/, which RABBIT_documentation.tex includes
% via \graphicspath{{figures/}}. Re-run this after the reference gait changes,
% then recompile the document (pdflatex twice, for the cross-references):
%
%   make_doc_figures
%   !cd docs && pdflatex RABBIT_documentation.tex
%
% Figures (all from the newest Results/col_result_*.mat):
%   fig_gait_stickfigure.png   one step, snapshots at N nodes      (Sec. collocation)
%   fig_col_vs_closedloop.png  collocation vs PD closed loop       (Sec. collocation)
%   fig_limit_cycle.png        phase portrait over valid steps     (Sec. stability)
%   fig_joint_torque.png       node angles + torques vs the box    (Sec. physical)
%
% The stick-figure and torque plots read the collocation NODE values directly
% (no simulation); the crosscheck and phase portrait run closed-loop sims and
% therefore take ~30 s in total.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

figdir = fullfile(root,'docs','figures');
if ~exist(figdir,'dir'), mkdir(figdir); end

% ---- newest collocation gait (name-sorted; see rabbit_hzd_collocation) ----
d = dir(fullfile(root,'Results','col_result_*.mat'));
assert(~isempty(d), 'No Results/col_result_*.mat');
[~,ord] = sort({d.name});
colfile = fullfile(d(ord(end)).folder, d(ord(end)).name);
fprintf('collocation gait: %s\n', d(ord(end)).name);

S = load(colfile);
p = S.p;
[X, U, ~, coeffs, T] = col_unpack(S.z_opt, p);
nq = p.nq;  N = p.N;
tnode = linspace(0, T, N);

names = {'q_1 stance hip','q_2 stance knee','q_3 swing hip','q_4 swing knee'};

%% ---------------- Fig 1: collocation vs closed loop --------------------
% col_crosscheck already draws exactly this overlay; reuse it so the figure
% cannot drift from the diagnostic.
col_crosscheck(colfile, 'pd', true);
exportgraphics(gcf, fullfile(figdir,'fig_col_vs_closedloop.png'), 'Resolution',150);
close(gcf);
fprintf('  [1/4] fig_col_vs_closedloop.png\n');

%% ---------------- Fig 2: stick figure sequence -------------------------
f = figure('Color','w','Position',[100 100 900 320]);
hold on; grid on; axis equal;
nsnap = 7;
idx = round(linspace(1, N, nsnap));
shift = 0;                       % lay snapshots out left-to-right
dx = 0.16;
for k = 1:nsnap
    q = X(1:nq, idx(k));
    pts = get_body_points(q);
    sf = pts.stance_foot; wf = pts.swing_foot; hp = pts.hip;
    sk = pts.stance_knee; wk = pts.swing_knee; tt = pts.torso_top;
    o = shift;
    a = 0.25 + 0.75*(k-1)/(nsnap-1);          % fade in over the step
    plot([hp(1) sk(1) sf(1)]+o, [hp(2) sk(2) sf(2)], '-', ...
         'Color',[0 0.2 0.8 a], 'LineWidth',3);
    plot([hp(1) wk(1) wf(1)]+o, [hp(2) wk(2) wf(2)], '-', ...
         'Color',[0.85 0.1 0.1 a], 'LineWidth',3);
    plot([hp(1) tt(1)]+o, [hp(2) tt(2)], '-', ...
         'Color',[0 0.55 0.2 a], 'LineWidth',4);
    plot(hp(1)+o, hp(2), 'ko','MarkerFaceColor','k','MarkerSize',5);
    shift = shift + dx;
end
plot([-0.3 shift+0.5],[0 0],'k-','LineWidth',1.5);
xlabel('x [m]  (snapshots offset for clarity)'); ylabel('z [m]');
title(sprintf('One collocation step, %d snapshots  (blue = stance, red = swing, green = torso)', nsnap));
ylim([-0.05 1.85]);
exportgraphics(f, fullfile(figdir,'fig_gait_stickfigure.png'), 'Resolution',150);
close(f);
fprintf('  [2/4] fig_gait_stickfigure.png\n');

%% ---------------- Fig 3: joint angles + torques ------------------------
umax = getfielddef(p,'u_max',120);
enf  = getfielddef(p,'enforce_torque',false);
f = figure('Color','w','Position',[100 100 860 620]);

subplot(2,1,1); hold on; grid on;
co = lines(4);
for j = 1:4
    plot(tnode, X(3+j,:), '-', 'Color',co(j,:), 'LineWidth',1.6, ...
         'DisplayName',names{j});
end
xlabel('t [s]'); ylabel('angle [rad]');
title('Actuated joint trajectories at the collocation nodes');
legend('Location','northwest','NumColumns',2,'FontSize',8);

subplot(2,1,2); hold on; grid on;
for j = 1:4
    plot(tnode, U(j,:), '-', 'Color',co(j,:), 'LineWidth',1.6, ...
         'DisplayName',names{j});
end
yline( umax, 'k--', 'LineWidth',1.3, 'DisplayName',sprintf('\\pm u_{max} = %g Nm',umax));
yline(-umax, 'k--', 'LineWidth',1.3, 'HandleVisibility','off');
% Mark every node torque that leaves the box -- with enforce_torque off the
% solver is free to exceed it, and this is where it does.
[jr, jc] = find(abs(U) > umax);
if ~isempty(jr)
    plot(tnode(jc), U(sub2ind(size(U),jr,jc)), 'ro', 'MarkerSize',11, ...
         'LineWidth',1.8, 'DisplayName','violates box');
end
xlabel('t [s]'); ylabel('torque [Nm]');
title(sprintf('Node torques vs the Table 3.1 box  (peak |u| = %.1f Nm, enforce\\_torque = %d)', ...
      max(abs(U(:))), enf));
ylim([-1.45*umax 1.15*umax]);
legend('Location','southeast','NumColumns',3,'FontSize',8);
exportgraphics(f, fullfile(figdir,'fig_joint_torque.png'), 'Resolution',150);
close(f);
fprintf('  [3/4] fig_joint_torque.png  (peak |u| = %.2f Nm, limit %g)\n', max(abs(U(:))), umax);

%% ---------------- Fig 4: limit cycle over several steps ----------------
% The seed orbit is orbitally UNSTABLE (rho ~ 3), so the closed loop departs
% within a few steps: from step 3 on, |q| reaches ~700 rad, which is numerical
% blow-up rather than walking. Plot only the steps that stay physically
% meaningful -- the interesting content is that the orbit does not CLOSE.
QLIM = 5;                                   % rad; beyond this it is blow-up
p.controller = 'pd';
x_cur = X(:,1);
segs = {};  starts = [];
for s = 1:6
    try
        [~, xs] = simulate_hzd_gait_full(coeffs, x_cur, p);
    catch
        break;
    end
    if isempty(xs) || ~all(isfinite(xs(:))), break; end
    if max(abs(xs(:,1:nq)),[],'all') > QLIM
        fprintf('  step %d exceeds %g rad (blow-up); not plotted\n', s, QLIM);
        break;
    end
    segs{end+1} = xs;  starts(:,end+1) = xs(1,:).'; %#ok<AGROW>
    x_cur = rabbit_reset_map(rabbit_impact_map(xs(end,:).'));
    if ~all(isfinite(x_cur)), break; end
end
nok = numel(segs);
fprintf('  limit-cycle: %d physically-valid steps plotted\n', nok);

f = figure('Color','w','Position',[100 100 880 680]);
cs = [0 0.35 0.75; 0.85 0.25 0.1; 0.2 0.6 0.2];
for j = 1:4
    subplot(2,2,j); hold on; grid on;
    for s = 1:nok
        xs = segs{s};
        plot(xs(:,3+j), xs(:,nq+3+j), '-', 'Color',cs(min(s,3),:), ...
             'LineWidth',1.3, 'DisplayName',sprintf('step %d',s));
    end
    % Poincare section: where each step BEGINS. If the gait were an orbitally
    % stable limit cycle these markers would coincide.
    for s = 1:nok
        plot(starts(3+j,s), starts(nq+3+j,s), 'o', 'MarkerSize',7, ...
             'MarkerFaceColor',cs(min(s,3),:), 'MarkerEdgeColor','k', ...
             'HandleVisibility','off');
    end
    xlabel(sprintf('%s [rad]', names{j})); ylabel('velocity [rad/s]');
    title(names{j});
    if j==1, legend('Location','best','FontSize',8); end
end
sgtitle(sprintf(['Phase portrait, %d PD closed-loop steps (circles = step starts).\n' ...
                 'The starts do not coincide: the orbit is not attracting.'], nok));
exportgraphics(f, fullfile(figdir,'fig_limit_cycle.png'), 'Resolution',150);
close(f);
fprintf('  [4/4] fig_limit_cycle.png\n');

fprintf('\nAll figures written to docs/figures/\n');
end

function v = getfielddef(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
