function fig = ch6_plot_barriers(sim, p, file)
%CH6_PLOT_BARRIERS  The barrier constraints along a run.  Figs. 6.13-6.15, 6.25(b).
%
%   fig = ch6_plot_barriers(sim, p)
%   fig = ch6_plot_barriers(sim, p, file)
%
% Three panels, and they are three DIFFERENT claims. The chapter's figures show
% only the first; the other two are what tell you whether the first is holding
% by design or by luck.
%
%   1. g_i(x) >= 0   -- the position constraint itself, (6.9)/(6.20). This is
%      Fig. 6.13-6.15 and Fig. 6.25(b). It is what the robot must satisfy.
%
%   2. h_CBF = gamma_b g + gdot >= 0   -- the relative-degree-1 function of
%      (6.1). This is the quantity the controller can actually act on, and it
%      goes negative BEFORE g does. A run where panel 1 is fine and panel 2 dips
%      is one where the guarantee lapsed and the constraint survived anyway --
%      worth knowing, and invisible in the chapter's own figures.
%
%   3. the QP row margin -- gddot + Kb eta_b evaluated at the u actually
%      applied. Zero means the row was tight and the barrier was doing work;
%      negative means the QP could not meet it, which under a hard row means it
%      was infeasible and fell back. Samples where that happened are marked.
%
% Panel 2 dipping while panel 1 holds is the normal state of affairs near
% touchdown and is not a failure: g -> 0 as the foot lands on the edge of a
% small stone while it is still moving at metres per second, so g/|gdot| is a
% few milliseconds there and no finite gamma_b keeps h_CBF positive at the very
% last instant. That is a property of the CONSTRUCTION, not of the tuning --
% see docs/CH6_STEPPING_STONES.md.
%
% Inputs
%   sim  : output of ch6_simulate
%   p    : parameter struct
%   file : optional PNG path
%
% See also CH6_PLOT_STONES, CH6_REPORT, CH6_CBF_ROW.

if nargin < 3, file = ''; end

n = sim.n_ok;
if n == 0 || isempty(sim.steps(1).log.h)
    error('ch6_plot_barriers:empty', ...
          'Nothing to plot: no completed steps, or the run carried no barriers.');
end

nb = size(sim.steps(1).log.h, 1);
B0 = ch6_barrier(sim.steps(1).x(:,1), [], with_stone(p, sim.steps(1)), 0);
labels = arrayfun(@(b) b.label, B0, 'UniformOutput', false);

cols = lines(max(nb, 3));

fig = figure('Color', 'w', 'Position', [80 60 980 860]);
% TeX, not LaTeX. MATLAB's default 'tex' interpreter handles \gamma and _{} but
% has no math mode, so a '$...$' fragment is passed through literally; asking
% for 'latex' instead makes the PLAIN parts of the string illegal (a bare '--',
% an unescaped '('), which is what the "valid interpreter syntax" warning was.
% Keeping every label in one interpreter is what avoids having to know which.
names = {'g_i(x)  --  the position constraint', ...
         'h_{CBF} = \gamma_b g + dg/dt  --  what the QP acts on', ...
         'row margin: d^2g/dt^2 + K_b \eta_b at the applied u'};

for panel = 1:3
    subplot(3,1,panel); hold on; box on;
    t_off = 0;
    for k = 1:n
        lg = sim.steps(k).log;
        tt = t_off + lg.t;
        switch panel
            case 1, Y = lg.h;
            case 2, Y = lg.h_cbf;
            case 3, Y = lg.margin;
        end
        for j = 1:nb
            plot(tt, Y(j,:), '-', 'Color', cols(j,:), 'LineWidth', 1.6);
        end
        if panel == 3 && any(~lg.feasible)
            plot(tt(~lg.feasible), zeros(1, sum(~lg.feasible)), 'kx', ...
                 'MarkerSize', 5, 'LineWidth', 1.2);
        end
        xline(t_off, ':', 'Color', [0.6 0.6 0.6]);
        t_off = t_off + sim.steps(k).T;
    end
    yline(0, 'k--', 'LineWidth', 1);
    xlabel('Time (s)');
    title(names{panel});
    if panel == 1
        legend(labels, 'Location', 'best', 'Interpreter', 'none');
    end
    if panel == 3, set(gca, 'YScale', 'linear'); end
    grid on;
end

sgtitle(sprintf('CBF constraints  --  %s barrier, poles (%.4g, %.4g)', ...
                p.cbf.form, p.cbf.gamma_b, p.cbf.gamma));

if ~isempty(file)
    exportgraphics(fig, file, 'Resolution', 150);
    fprintf('[ch6_plot_barriers] wrote %s\n', file);
end

end

function q = with_stone(p, step)
q = p;
q.stone = step.stone;
end
