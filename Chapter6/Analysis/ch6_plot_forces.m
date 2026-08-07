function fig = ch6_plot_forces(sim, p, file)
%CH6_PLOT_FORCES  Contact forces and torques along a run.  Fig. 6.26.
%
%   fig = ch6_plot_forces(sim, p)
%   fig = ch6_plot_forces(sim, p, file)
%
% The three constraints of (6.25) plus input saturation, exactly as Fig. 6.26
% draws them:
%
%   (a) normal force   Fv_st >= delta_N          -- keep the foot loaded
%   (b) friction cone  |Fh_st / Fv_st| <= kf     -- do not slip
%   (c) control inputs |u| <= u_max              -- do not saturate
%
% Every limit line is drawn whether or not it is enforced, and the line style
% says which: SOLID for a limit the QP is carrying as a row, DASHED for one
% that is only being measured. A dashed line that the data crosses is not a bug
% -- it is the measure-then-tighten workflow showing what enabling that row
% would cost.
%
% The forces come from lambda = lam_drift + lam_in u, the same KKT solve the
% controller used, so panel (a) and (b) are the forces the simulated robot
% applied rather than a reconstruction from the trajectory.
%
% Inputs
%   sim  : output of ch6_simulate
%   p    : parameter struct
%   file : optional PNG path
%
% See also CH6_REPORT, CH6_PLOT_BARRIERS, CH3_FORCES.

if nargin < 3, file = ''; end

n = sim.n_ok;
if n == 0
    error('ch6_plot_forces:empty', 'The run completed no steps.');
end

fig = figure('Color', 'w', 'Position', [80 60 980 860]);

panels = { 'F^v_{st}  (N)', '|F^h_{st} / F^v_{st}|', 'u  (Nm)' };
lims   = [p.limits.Fz_min, p.limits.mu_s, p.limits.u_max];
on     = [p.limits.enable.grf, p.limits.enable.friction, p.limits.enable.torque];

for panel = 1:3
    subplot(3,1,panel); hold on; box on;
    t_off = 0;
    for k = 1:n
        lg = sim.steps(k).log;
        tt = t_off + lg.t;
        switch panel
            case 1, plot(tt, lg.Fz, '-', 'Color', [0.10 0.35 0.75], 'LineWidth', 1.6);
            case 2
                m = lg.mu_fric;
                m(lg.Fz <= 1e-6) = NaN;    % |Fx/Fz| is meaningless at Fz = 0
                plot(tt, m, '-', 'Color', [0.75 0.25 0.10], 'LineWidth', 1.6);
            case 3, plot(tt, lg.u.', '-', 'LineWidth', 1.2);
        end
        xline(t_off, ':', 'Color', [0.6 0.6 0.6]);
        t_off = t_off + sim.steps(k).T;
    end

    sty = '--';  if on(panel), sty = '-'; end
    yline(lims(panel), sty, 'Color', [0.1 0.6 0.1], 'LineWidth', 1.4);
    if panel == 3
        yline(-lims(3), sty, 'Color', [0.1 0.6 0.1], 'LineWidth', 1.4);
    end

    xlabel('Time (s)');
    ylabel(panels{panel});
    title(sprintf('(%c)   limit %.4g   %s', 'a'+panel-1, lims(panel), ...
                  enf(on(panel))));
    grid on;
end

sgtitle('Ground reaction, friction cone and input saturation   (6.25)');

if ~isempty(file)
    exportgraphics(fig, file, 'Resolution', 150);
    fprintf('[ch6_plot_forces] wrote %s\n', file);
end

end

function s = enf(b)
if b, s = '[ENFORCED as a QP row]'; else, s = '(measured only)'; end
end
