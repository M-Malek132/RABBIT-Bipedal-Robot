function V = ch3_col_verify(z, p, verbose)
%CH3_COL_VERIFY  Is the collocation solution an ACTUAL trajectory?
%
%   V = ch3_col_verify(z, p)
%   V = ch3_col_verify(z, p, verbose)
%
% THIS IS THE CHECK THAT KEEPS DIRECT COLLOCATION HONEST, and it is not
% optional.
%
% Small defects mean the discrete equations are satisfied. They do NOT mean
% the node states approximate a solution of the ODE. When the mesh is too
% coarse for the dynamics, the optimizer can find a SPURIOUS DISCRETE SOLUTION
% -- defects at 1e-7 while the nodes sit far from any real trajectory -- and
% every number derived from it (step length, duration, speed, cost, the
% Table 3.1 quantities) is then fiction.
%
% Measured example from this project, at N = 15 nodes on a converged solve:
%
%       interval-1 defect                     7.18e-07
%       |x_2(collocation) - x_2(true flow)|   3.68e-02      <-- 5 orders worse
%
% The tell was that the accelerations reached 81 rad/s^2, so dq changed by
% ~2.3 within one h = 0.072 s interval. Hermite-Simpson's truncation error
% scales like h^5 times a fourth derivative; against dynamics that violent, 15
% nodes are simply not enough. The forward simulation duly disagreed with the
% collocation (T: 1.00 s vs 0.72 s), which is how the problem first surfaced.
%
% The check: integrate the closed-loop ODE from node 1 with a tight tolerance
% and compare against the node states at the same times. If the mesh is
% adequate the two agree to integration accuracy; if not, they diverge and the
% solution must be rejected and re-solved on a finer mesh.
%
% Inputs
%   z, p    : solution and parameters
%   verbose : print a per-node table (default false)
%
% Output
%   V : struct
%       .max_dev        worst |X_node - X_true| over all nodes
%       .dev            1 x N per-node deviation
%       .rel_dev        max_dev normalized by the state's own scale
%       .ok             true if max_dev <= p.verify_tol
%       .integrated     true if the reference rollout completed
%       .X_true         the reference trajectory at the node times
%
% See also CH3_COL_EVAL, CH3_REPORT, CH3_COL_SOLVE.

if nargin < 3, verbose = false; end

% Saved results carry the parameter struct that produced them, which may
% predate fields this function reads.
p = ch3_upgrade_params(p);

[X, T, alpha] = ch3_col_unpack(z, p);
N = size(X, 2);

% Pure feedforward, matching what the transcription assumes as its input.
pf = p;
pf.controller = 'ff';

opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);

V = struct('max_dev', Inf, 'dev', nan(1,N), 'rel_dev', Inf, ...
           'ok', false, 'integrated', false, 'X_true', []);

try
    sol = ode45(@(t,x) ch3_ode_rhs(t,x,alpha,pf), [0 T], X(:,1), opts);
catch
    if verbose
        fprintf('  ch3_col_verify: reference integration FAILED to start.\n');
    end
    return;
end

% ode45 can stop short if the dynamics become too stiff to step through --
% itself a signal that the mesh cannot be resolving them either.
t_avail = sol.x(end);
t_nodes = linspace(0, T, N);
usable  = t_nodes <= t_avail + 1e-12;

X_true = nan(size(X));
X_true(:, usable) = deval(sol, min(t_nodes(usable), t_avail));

dev = nan(1, N);
dev(usable) = max(abs(X(:, usable) - X_true(:, usable)), [], 1);

V.X_true    = X_true;
V.dev       = dev;
V.max_dev   = max(dev(usable));
V.integrated = all(usable);
scale       = max(1, max(abs(X(:))));
V.rel_dev   = V.max_dev / scale;

tol = p.verify_tol;
V.ok = V.integrated && isfinite(V.max_dev) && V.max_dev <= tol;

if verbose
    fprintf('\n COLLOCATION VERIFICATION (nodes vs a tight rollout from node 1)\n');
    if ~V.integrated
        fprintf('   reference rollout stopped early at t = %.4f of T = %.4f\n', t_avail, T);
        fprintf('   -- the dynamics are too stiff to integrate; the mesh cannot be\n');
        fprintf('      resolving them either. REJECT this solution.\n');
    end
    fprintf('   max |X_node - X_true|        %.3e   (tol %.1e)  -> %s\n', ...
            V.max_dev, tol, ternary(V.ok, 'OK', 'REJECT: mesh too coarse'));
    if ~V.ok
        fprintf('   The defects can be tiny while this is large: that is a SPURIOUS\n');
        fprintf('   discrete solution. Re-solve with more nodes (p.N_nodes).\n');
    end
    k_show = unique(round(linspace(1, N, min(N, 8))));
    fprintf('   node : ');   fprintf('%9d', k_show);   fprintf('\n');
    fprintf('   dev  : ');   fprintf('%9.1e', dev(k_show)); fprintf('\n');
end

end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end
