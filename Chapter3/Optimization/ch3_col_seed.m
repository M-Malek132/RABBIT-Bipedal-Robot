function [z0, info] = ch3_col_seed(p, x0)
%CH3_COL_SEED  Build an initial collocation guess by rolling out the seed gait.
%
%   [z0, info] = ch3_col_seed(p)
%   [z0, info] = ch3_col_seed(p, x0)
%
% Builds z0 = [X; T; alpha] by integrating the seed under the FEEDFORWARD
% controller and sampling nodes from the result.
%
% Why simulate rather than interpolate poses:  ch3_seed puts x0 exactly on the
% zero dynamics surface (y = 0, ydot = 0), and u_ff renders that surface
% invariant, so the rolled-out trajectory stays on Z for the whole step. The
% resulting guess therefore satisfies the node-1 equalities AND the defects to
% integration accuracy from the very first fmincon evaluation. The only
% conditions left genuinely violated are the ones a seed cannot know --
% periodicity, and the NEC1 walking rate. That is a far better starting point
% than interpolating between poses, which violates everything at once.
%
% NODES ARE SAMPLED ON THE PHASE, NOT ON TIME.  The step is cut where theta
% first reaches theta_plus, so theta(x_N) = theta_plus holds at the seed too.
% Sampling a fixed duration instead would leave the terminal phase equality
% violated by an arbitrary amount.
%
% Inputs
%   p  : parameter struct
%   x0 : optional start state (default: ch3_seed's hand-tuned pose)
%
% Outputs
%   z0   : packed decision vector
%   info : struct .x0 .alpha .T .theta_reached .guard_fired .sol
%
% See also CH3_SEED, CH3_COL_SOLVE.

if nargin < 2, x0 = []; end

[alpha, x0] = ch3_seed(p, x0);

% Roll out under pure feedforward: no feedback term to perturb the zero
% dynamics surface the seed was constructed to lie on.
p_ff = p;
p_ff.controller = 'ff';

s = ch3_step(x0, alpha, p_ff);

% Find where theta first reaches theta_plus.
th = p.c_theta(:).' * s.x(1:p.nq, :);
idx = find(th >= p.theta_plus, 1, 'first');

theta_reached = th(end);
if isempty(idx)
    % The seed never swept the full phase. Use the whole step and let fmincon
    % pull theta_end onto theta_plus -- reported so the caller knows the seed
    % starts further from feasible than usual.
    t_end = s.t(end);
else
    % Linear interpolation onto the exact crossing.
    if idx == 1
        t_end = s.t(1);
    else
        w = (p.theta_plus - th(idx-1)) / (th(idx) - th(idx-1));
        t_end = s.t(idx-1) + w*(s.t(idx) - s.t(idx-1));
    end
    theta_reached = p.theta_plus;
end

t_end = max(t_end, p.T_min);
T     = min(max(t_end, p.T_min), p.T_max);

N     = p.N_nodes;
t_nodes = linspace(0, T, N);

% Resample with DEVAL, not interp1.  deval uses ode45's own dense-output
% polynomial, so the node states are accurate to the integration tolerance
% (~1e-9 here).  Re-interpolating the returned grid with pchip instead leaves
% an interpolation error around 4e-6 baked into the nodes, which shows up as a
% defect floor that does not shrink when N grows -- the transcription then
% looks like it has lost its order of accuracy when in fact the seed was the
% inaccurate part.
X = deval(s.sol, t_nodes);

z0 = ch3_col_pack(X, T, alpha, p);

info = struct('x0', x0, 'alpha', alpha, 'T', T, ...
              'theta_reached', theta_reached, ...
              'guard_fired', s.ok, 'sol', s);

end
