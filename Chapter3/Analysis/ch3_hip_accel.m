function [fig, A] = ch3_hip_accel(z, p, n_steps, savepath)
%CH3_HIP_ACCEL  Hip acceleration along the periodic motion of a gait.
%
%   [fig, A] = ch3_hip_accel(z, p)
%   [fig, A] = ch3_hip_accel(z, p, n_steps, 'Results/ch3_hip_accel.png')
%
% Simulates n_steps of the hybrid system from node 1 of the collocation
% solution and reports the world-frame acceleration of the hip.
%
% THE ACCELERATION IS COMPUTED FROM THE DYNAMICS, NOT DIFFERENCED IN TIME.
%
%       a = J(q) ddq + (dJ/dt) dq,    ddq = the bottom half of f(x) + g(x) u(x)
%
% That distinction is not cosmetic here.  The hybrid motion has a velocity
% JUMP at every foot strike -- the impact is impulsive, so the acceleration
% there is a Dirac, not a large finite number.  Any finite difference taken
% across a strike reports a spike whose height is set by the sample spacing
% rather than by the gait, and refining the grid makes it taller.  Computing
% a pointwise from the closed loop sidesteps that entirely; each step is
% drawn as its own curve, so the impulsive part shows up as the break it is.
%
% (dJ/dt)dq is obtained as the directional derivative d/de[J(q + e dq) dq] at
% e = 0 -- one more central difference, rather than assembling the full
% second derivative of the forward kinematics.
%
% Panel 4 folds every step onto tau = t/T_step.  Overlapping curves are a
% direct read of periodicity in the quantity being plotted, which is stricter
% than the periodicity residual: ceq constrains the STATE at the step
% boundary, while a depends on ddq and so on the control as well.
%
% Inputs
%   z        : collocation decision vector
%   p        : parameter struct
%   n_steps  : steps to simulate (default 6)
%   savepath : optional PNG path
%
% Outputs
%   fig : figure handle
%   A   : struct with .t .tau .a .v .pos per step, plus
%           .spread   max spread of a across steps on the tau grid [m/s^2]
%           .fd_check max |a - dv/dt| inside a step, the validation residual
%           .fd_rel   the same, relative to the peak |a| on that step
%           .peak     max |a| [m/s^2]
%
% See also CH3_PLOT_GAIT, CH3_FORCES, CH3_BODY_POINTS.

if nargin < 3 || isempty(n_steps),  n_steps = 6;   end
if nargin < 4,                      savepath = ''; end

NPT = 400;                          % samples per step

[X, ~, alpha] = ch3_col_unpack(z, p);
sim = ch3_simulate(X(:,1), alpha, p, n_steps);

if sim.n_ok < 1
    error('ch3_hip_accel:nosteps', 'No step completed: %s', sim.reason);
end
if sim.failed
    warning('ch3_hip_accel:partial', '%d/%d steps only: %s', ...
            sim.n_ok, n_steps, sim.reason);
end

%% per-step kinematics on a uniform grid inside each step
NS = sim.n_ok;
St = cell(1, NS);
t_off = 0;
for k = 1:NS
    s  = sim.steps(k);
    tt = linspace(0, s.T, NPT);

    % deval resamples at the solver's own accuracy; interp1 on the output
    % grid would be a second, coarser approximation on top of the first.
    if ~isempty(s.sol)
        XX = deval(s.sol, tt);
    else
        XX = interp1(s.t.', s.x.', tt.').';      % sampled-data branch
    end

    a = zeros(2, NPT); v = zeros(2, NPT); pos = zeros(2, NPT);
    for j = 1:NPT
        x   = XX(:,j);
        q   = x(1:p.nq);
        dq  = x(p.nq+1:end);
        xd  = ch3_ode_rhs(0, x, alpha, p);
        [pos(:,j), v(:,j), a(:,j)] = hip_kin(q, dq, xd(p.nq+1:end));
    end

    St{k} = struct('t', t_off + tt, 'tau', tt/s.T, 'a', a, 'v', v, ...
                   'pos', pos, 'T', s.T, 'L', s.L_step);
    t_off = t_off + s.T;
end

%% validation: a against a time difference of v, away from the strikes
s1   = St{min(2, NS)};
dt   = s1.t(2) - s1.t(1);
a_fd = (s1.v(:,3:end) - s1.v(:,1:end-2)) / (2*dt);
fd_check = max(abs(a_fd - s1.a(:,2:end-1)), [], 'all');
fd_rel   = fd_check / max(vecnorm(s1.a));

%% periodicity of the acceleration itself, on the shared tau grid
allA = cell2mat(cellfun(@(s) s.a, St, 'UniformOutput', false));
if NS >= 2
    P = cell2mat(cellfun(@(s) s.a(:), St, 'UniformOutput', false));
    spread = max(max(P,[],2) - min(P,[],2));
else
    spread = NaN;
end

A = struct('steps', {St}, 'spread', spread, 'fd_check', fd_check, ...
           'fd_rel', fd_rel, 'peak', max(vecnorm(allA)), ...
           'T', cellfun(@(s) s.T, St), 'L', cellfun(@(s) s.L, St));

fprintf('hip acceleration: %d steps, peak |a| = %.1f m/s^2 (%.1f g)\n', ...
        NS, A.peak, A.peak/p.g0);
fprintf('  a vs dv/dt inside a step : %.2e m/s^2  (%.1e relative)\n', fd_check, fd_rel);
fprintf('  spread across steps      : %.2e m/s^2\n', spread);

%% figure
strikes = cumsum(A.T);
fig = figure('Color','w','Position',[60 60 1150 1050]);
co  = lines(7);
pad = @(y) [min(y) - 0.08*range(y), max(y) + 0.08*range(y)];

names = {'a_x  [m/s^2]', 'a_z  [m/s^2]', '|a|  [m/s^2]'};
ttl   = {'Hip horizontal acceleration', ...
         'Hip vertical acceleration   (+ up; the dashed line is free fall)', ...
         'Hip acceleration magnitude'};
for i = 1:3
    subplot(4,1,i); hold on; grid on
    for k = 1:NS
        if i < 3, yk = St{k}.a(i,:); else, yk = vecnorm(St{k}.a); end
        plot(St{k}.t, yk, '-', 'Color', co(i,:), 'LineWidth', 1.5);
    end
    for k = 1:NS-1, xl = xline(strikes(k)); xl.Color = [.75 .75 .75]; end
    yline(0, 'k-');
    switch i
        case 1
            ylim(pad(allA(1,:)));
        case 2
            yl = yline(-p.g0, '--'); yl.Color = [.6 0 .6];
            text(0.015, -p.g0, ' -g', 'Color', [.6 0 .6], ...
                 'VerticalAlignment', 'bottom', 'FontSize', 8);
            ylim(pad([allA(2,:), -1.07*p.g0]));
        case 3
            yl = yline(p.g0, ':'); yl.Color = [.6 0 .6];
            text(0.015, p.g0, ' g', 'Color', [.6 0 .6], ...
                 'VerticalAlignment', 'bottom', 'FontSize', 8);
            ylim(pad(vecnorm(allA)));
            xlabel('t  [s]');
    end
    ylabel(names{i}); title(ttl{i}, 'FontWeight', 'bold');
    xlim([0 strikes(end)]);
end

subplot(4,1,4); hold on; grid on
for k = 1:NS
    plot(St{k}.tau, St{k}.a(1,:), '-', 'Color', [co(1,:) 0.7], 'LineWidth', 1.2);
    plot(St{k}.tau, St{k}.a(2,:), '-', 'Color', [co(2,:) 0.7], 'LineWidth', 1.2);
end
xlabel('\tau = t / T_{step}'); ylabel('a  [m/s^2]');
title(sprintf(['All %d steps folded onto one period  ' ...
               '(they coincide to %.1e m/s^2 -- the motion is periodic)'], NS, spread), ...
      'FontWeight', 'bold');
legend({'a_x','a_z'}, 'Location', 'northwest', 'FontSize', 8);
xlim([0 1]); ylim(pad(allA(:).'));

sgtitle(sprintf(['Hip acceleration over the periodic motion  |  %s, %s basis\n' ...
                 '%d steps  |  T = %.4f s, L = %.4f m, v = %.3f m/s  |  peak |a| = %.1f m/s^2 (%.1f g)\n' ...
                 'grey = foot strikes; velocity jumps there, so the impulsive part is a break, not a spike'], ...
        p.controller, p.basis, NS, mean(A.T), mean(A.L), mean(A.L)/mean(A.T), ...
        A.peak, A.peak/9.81), 'FontSize', 10);

if ~isempty(savepath)
    exportgraphics(fig, savepath, 'Resolution', 150);
    fprintf('Wrote %s\n', savepath);
end

end

% ---------------------------------------------------------------------------
function [pos, v, a] = hip_kin(q, dq, ddq)
%HIP_KIN  World-frame hip position, velocity and acceleration.
h  = 1e-6;
nq = numel(q);

pos = hip_pos(q);

J = zeros(2, nq);
for i = 1:nq
    qp = q; qp(i) = qp(i) + h;
    qm = q; qm(i) = qm(i) - h;
    J(:,i) = (hip_pos(qp) - hip_pos(qm)) / (2*h);
end
v = J * dq;

% (dJ/dt)dq = d/de [ J(q + e dq) dq ] at e = 0. The step is scaled by |dq|
% so the perturbation stays the same size in q whatever the velocity is.
e = h / max(1, norm(dq));
a = J * ddq + (jac_along(q + e*dq, dq, h) - jac_along(q - e*dq, dq, h)) / (2*e);
end

% ---------------------------------------------------------------------------
function w = jac_along(q, dq, h)
%JAC_ALONG  J(q) dq, without forming J.
w = zeros(2,1);
for i = 1:numel(q)
    qp = q; qp(i) = qp(i) + h;
    qm = q; qm(i) = qm(i) - h;
    w = w + (hip_pos(qp) - hip_pos(qm)) / (2*h) * dq(i);
end
end

% ---------------------------------------------------------------------------
function pos = hip_pos(q)
b   = ch3_body_points(q);
pos = b.hip;
end
