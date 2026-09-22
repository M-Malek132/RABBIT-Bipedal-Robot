function [z, out] = ch6_lib_newton(p, L_target, z0, opts)
%CH6_LIB_NEWTON  One library gait by Newton continuation on Chapter 3's constraints.
%
%   [z, out] = ch6_lib_newton(p, L_target, z0)
%   [z, out] = ch6_lib_newton(p, L_target, z0, opts)
%
% The corrector of a predictor-corrector march in step length. It solves
%
%       ceq(z) = 0                   Chapter 3's defects, periodicity, outputs
%       L_step(z) - L_target = 0     the one row Section 6.4 adds
%       c_i(z) = 0,  i active        inequalities that are binding or violated
%
% by Levenberg-Marquardt steps (below), keeping z inside
% ch3_col_bounds. Nothing about the transcription changes -- the rows are
% ch3_col_constraints' own -- so a converged z is a Chapter-3 gait and is put
% through ch3_col_verify like one.
%
% WHY NOT FMINCON, WHICH CH6_LIB_SOLVE USES. Measured on posture_195's 61-node
% mesh (879 variables, 869 equalities): with Chapter 3's torque cost SQP left
% the feasible region in two iterations (the seed is not a minimiser of that
% cost, optimality 8.4e4), and with a proximal cost it held feasibility but took
% line-search steps of 2e-4 to 6e-2 at about two minutes per iteration -- hours
% per 1 cm rung. A library needs neighbouring FEASIBLE gaits, not re-optimised
% ones, and Newton's method on the constraints converges in a handful of
% Jacobians from a 1 cm predictor. The minimum-norm step is also the smallest
% change to z that fixes the residual, which is what (6.22) wants between
% neighbours.
%
% THE JACOBIAN is forward differences, split into p.lib.workers chunks on the
% current pool (33 s for 879 columns on 6 workers here), serial without one.
%
% Inputs
%   p        : the SEED's parameter struct (ch6_lib_build passes it)
%   L_target : step length [m]
%   z0       : predictor
%   opts     : .tol (1e-6) .max_iter (25) .h (1e-6) .act_tol (1e-4)
%              .ineq (false) carry active inequalities as rows
%
% Outputs
%   z   : corrected decision vector
%   out : .ok .iters .res .L_step .T .speed .alpha .J_ch3 .verify .wall_time
%
% See also CH6_LIB_BUILD, CH6_LIB_SOLVE, CH3_COL_CONSTRAINTS, CH3_COL_VERIFY.

if nargin < 4, opts = struct(); end
d = struct('tol', 1e-6, 'max_iter', 25, 'h', 1e-6, 'act_tol', 1e-4, 'ineq', false);
for f = fieldnames(d).'
    if ~isfield(opts, f{1}), opts.(f{1}) = d.(f{1}); end
end

p = ch3_upgrade_params(p);
p.enforce_nec1 = false;

N = (numel(z0) - 1 - p.ny*p.n_ctrl) / p.nx;
if isfield(p, 'lib') && isfield(p.lib, 'T_band') && ~isempty(p.lib.T_band)
    [~, T0] = ch3_col_unpack(z0, p);
    p.T_min = max(p.T_min, T0 * (1 - p.lib.T_band));
    p.T_max = min(p.T_max, T0 * (1 + p.lib.T_band));
end
[lb, ub] = ch3_col_bounds(p, N);

t0 = tic;
z  = min(max(z0, lb), ub);
[r, act] = residual(z, p, L_target, opts.act_tol, [], opts.ineq);
fprintf('[ch6_lib_newton] L* = %.4f  it 0  |r| %.2e  (%d active ineq)\n', ...
        L_target, max(abs(r)), nnz(act));

% LEVENBERG-MARQUARDT, NOT PLAIN NEWTON. Measured at the seed: the equality
% Jacobian is rank deficient (J J' singular -- the rows are not independent),
% so the minimum-norm Newton step is |dz| = 342, and a few rows depend so
% nonlinearly on the Bezier coefficients (the feedforward goes through
% LgLf y, which is poorly conditioned at some nodes) that the residual along
% that step grows at every step length tried down to 1e-3. The damped step
% dz = -(J'J + mu I)^-1 J' r is well defined whatever the rank and shrinks
% toward steepest descent as mu grows; a rejected trial costs one residual
% evaluation, not a Jacobian.
mu = 1e-2;
it = 0;  n_trial = 0;
while max(abs(r)) > opts.tol && it < opts.max_iter
    it = it + 1;
    Jf = jacobian(z, p, L_target, opts.h);          % all rows, c then ceq then L
    J  = Jf(row_mask(act, size(Jf, 1)), :);
    JtJ = J.' * J;  Jtr = J.' * r;
    scale = max(diag(JtJ));

    base = norm(r);  improved = false;
    while mu <= 1e8
        dz = -(JtJ + mu * scale * eye(numel(z))) \ Jtr;
        zt = min(max(z + dz, lb), ub);
        rt = residual(zt, p, L_target, opts.act_tol, act, opts.ineq);
        n_trial = n_trial + 1;
        if all(isfinite(rt)) && norm(rt) < base
            improved = true;  mu = max(mu / 3, 1e-12);
            break;
        end
        mu = mu * 5;
    end
    if ~improved
        fprintf('[ch6_lib_newton] L* = %.4f  it %d  no decrease at any damping; stopping\n', ...
                L_target, it);
        break;
    end
    z = zt;
    [r, act] = residual(z, p, L_target, opts.act_tol, [], opts.ineq);
    fprintf('[ch6_lib_newton] L* = %.4f  it %d  mu %.1e  |dz| %.2e  |r| %.2e  (%d active)  %.0f s\n', ...
            L_target, it, mu, norm(dz), max(abs(r)), nnz(act), toc(t0));
end

[~, T, alpha] = ch3_col_unpack(z, p);
E = ch3_col_eval(z, p);
[c, ceq] = ch3_col_constraints(z, p);
out = struct('iters', it, 'res', max(abs(r)), ...
             'max_ceq', max(abs([ceq; E.L_step - L_target])), 'max_c', max(c), ...
             'L_step', E.L_step, 'T', T, 'speed', E.L_step / T, 'alpha', alpha, ...
             'L_target', L_target, 'J_ch3', ch3_col_cost(z, p), ...
             'fval', ch3_col_cost(z, p), 'exitflag', double(max(abs(r)) <= opts.tol), ...
             'wall_time', toc(t0), 'verify', ch3_col_verify(z, p, false));
out.ok = out.exitflag == 1 && out.max_c <= opts.tol && out.verify.ok;
fprintf(['[ch6_lib_newton] L* = %.4f -> L = %.4f, T = %.3f s, v = %.3f m/s, ' ...
         'ceq %.1e, c %.1e, mesh dev %.2e (%s), %.0f s\n'], L_target, out.L_step, ...
        T, out.speed, out.max_ceq, out.max_c, out.verify.max_dev, ...
        tf(out.verify.ok), out.wall_time);
end

% ---------------------------------------------------------------------------
function [r, act] = residual(z, p, L_target, act_tol, act, use_ineq)
[c, ceq] = ch3_col_constraints(z, p);
E = ch3_col_eval(z, p);
if ~use_ineq
    act = false(numel(c), 1);
elseif isempty(act)
    act = c(:) > -act_tol;
end
% An active inequality is driven to its bound; a satisfied one that the
% active set carried along contributes nothing once it is back inside.
ca = c(:);  ca = ca(act);
r  = [ca; ceq(:); E.L_step - L_target];
end

function m = row_mask(act, n_rows)
m = [act(:); true(n_rows - numel(act), 1)];
end

function J = jacobian(z, p, L_target, h)
f0 = all_rows(z, p, L_target);
n  = numel(z);
pool = gcp('nocreate');
if isempty(pool)
    J = zeros(numel(f0), n);
    for j = 1:n
        zz = z;  zz(j) = zz(j) + h;
        J(:, j) = (all_rows(zz, p, L_target) - f0) / h;
    end
    return;
end
nw = pool.NumWorkers;
cols = cell(1, nw);  blocks = cell(1, nw);
for w = 1:nw, cols{w} = w:nw:n; end
parfor w = 1:nw
    my = cols{w};
    B  = zeros(numel(f0), numel(my));
    for j = 1:numel(my)
        zz = z;  zz(my(j)) = zz(my(j)) + h;
        B(:, j) = (all_rows(zz, p, L_target) - f0) / h;
    end
    blocks{w} = B;
end
J = zeros(numel(f0), n);
for w = 1:nw, J(:, cols{w}) = blocks{w}; end
end

function f = all_rows(z, p, L_target)
[c, ceq] = ch3_col_constraints(z, p);
E = ch3_col_eval(z, p);
f = [c(:); ceq(:); E.L_step - L_target];
end

function s = tf(b)
if b, s = 'ok'; else, s = 'FAILED'; end
end
