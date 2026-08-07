function fig = ch6_phase_plot(sim, p, alpha_nom, file)
%CH6_PHASE_PLOT  Torso velocity against torso angle over a run.  Fig. 6.7.
%
%   fig = ch6_phase_plot(sim, p)
%   fig = ch6_phase_plot(sim, p, alpha_nom, file)
%
% Fig. 6.7 plots RABBIT's torso velocity against its torso angle over 30 random
% footholds, "the thick red line depicts the nominal limit cycle of the periodic
% walking gait for comparison". The point of the figure is that the stepping-
% stone trajectories stay in a bounded neighbourhood of the nominal orbit
% instead of wandering -- which is a stability statement the step-length plot
% cannot make.
%
% The torso is the UNACTUATED coordinate (p.iact = 4:7 excludes q_t), so its
% phase portrait is the zero-dynamics view: the four actuated joints are being
% driven wherever the barrier needs them, and this shows what the one degree of
% underactuation did in response. A plot of an actuated joint would mostly show
% the controller succeeding at tracking, which is not the question.
%
% The nominal cycle is drawn by simulating the SAME gait with no barrier and no
% terrain, so the comparison is against this repository's own limit cycle rather
% than against an idealisation of it.
%
% Inputs
%   sim       : output of ch6_simulate
%   p         : parameter struct
%   alpha_nom : the nominal gait; pass [] to skip the reference orbit
%   file      : optional PNG path
%
% See also CH6_SIMULATE, CH3_POINCARE, CH6_PLOT_STONES.

if nargin < 3, alpha_nom = []; end
if nargin < 4, file = ''; end

if sim.n_ok == 0
    error('ch6_phase_plot:empty', 'The run completed no steps.');
end

fig = figure('Color', 'w', 'Position', [100 100 760 620]);
hold on; box on;

%% ---------------------------------------------------- the nominal limit cycle
if ~isempty(alpha_nom)
    q = p;
    q.controller  = 'iolin_pd';
    q.cbf.problem = 'none';
    q.n_steps     = 6;
    x0 = sim.steps(1).x(:, 1);
    ref = ch6_simulate(x0, alpha_nom, q);
    if ref.n_ok >= 2
        % Skip the first step: it starts wherever the stepping-stone run
        % started, and the orbit is what the gait settles ONTO.
        k0 = ref.steps(2).t;  %#ok<NASGU>
        idx = size(ref.steps(1).x, 2) + 1;
        plot(ref.x(3, idx:end), ref.x(3+p.nq, idx:end), '-', ...
             'Color', [0.80 0.10 0.10], 'LineWidth', 3.5);
    end
end

%% ----------------------------------------------------- the stepping-stone run
plot(sim.x(3, :), sim.x(3+p.nq, :), '-', ...
     'Color', [0.15 0.30 0.70], 'LineWidth', 1.0);

% Impacts, where the velocity jumps.
t_off = 0;
for k = 1:sim.n_ok
    xe = sim.steps(k).x_end;
    plot(xe(3), xe(3+p.nq), 'k.', 'MarkerSize', 12);
    t_off = t_off + sim.steps(k).T;
end

xlabel('torso angle  q_t  (rad)');
ylabel('torso velocity  dq_t/dt  (rad/s)');
title(sprintf('Phase portrait over %d discrete footholds', sim.n_ok));
if ~isempty(alpha_nom)
    legend({'nominal limit cycle', 'stepping stones', 'impacts'}, ...
           'Location', 'best');
end
grid on;

if ~isempty(file)
    exportgraphics(fig, file, 'Resolution', 150);
    fprintf('[ch6_phase_plot] wrote %s\n', file);
end

end
