function ch5_gen_pendulum(varargin)
%CH5_GEN_PENDULUM  Derive and write the pendulum's exact Lie derivatives.
%
%   ch5_gen_pendulum()
%   ch5_gen_pendulum('optimize', false)   faster to generate, slower to run
%
% Writes Chapter5/Model/Generated/:
%
%   ch5_pend_theta_lie.m   theta^(2), theta^(3), and the (5.21) split of
%                          theta^(4) into  L_f^4 theta  +  L_g L_f^3 theta * tau
%   ch5_pend_py_lie.m      the same for the end-effector height py2, i.e. the
%                          eta_b stack of (5.10) plus the two rb-th order terms
%
% RUN ONCE, COMMIT THE OUTPUT. This mirrors what the repository already does
% for the RABBIT dynamics in Dynamics/: the Symbolic Math Toolbox is needed to
% REGENERATE, never to run. Regenerate only if ch5_pendulum or the parameter
% packing in ch5_pend_pv changes.
%
% ---------------------------------------------------- why generate at all
% The ECBF condition (5.9) needs L_f^rb h and L_g L_f^(rb-1) h EXACTLY. For
% rb = 4 on a 2-link arm those are fourth derivatives of a trigonometric
% expression through an inverted mass matrix, and the obvious alternatives are
% both bad:
%
%   * NESTED FINITE DIFFERENCES lose roughly three digits per level. By level 4
%     an h^(4) of order 1 carries an error of order 1e-3, and the ECBF row is
%     an INEQUALITY THAT IS MEANT TO BE TIGHT -- the controller sits exactly on
%     it whenever safety is active. An error there is not noise on a plot, it
%     is the constraint being enforced at the wrong level.
%   * COMPLEX-STEP is exact but does not nest: taking Im() at one level
%     destroys the perturbation the next level needs.
%
% Symbolic differentiation has neither problem, and moving it to generation
% time means the ODE right-hand side pays nothing for it.
%
% The generated code is checked against the live model in ch5_test_barrier, by
% differentiating along an actual drift trajectory with a high-order stencil --
% an independent route to the same numbers.
%
% See also CH5_PENDULUM, CH5_PEND_PV, CH5_BARRIER, CH5_TEST_BARRIER.

o = struct('optimize', true);
for i = 1:2:numel(varargin), o.(lower(varargin{i})) = varargin{i+1}; end

here    = fileparts(mfilename('fullpath'));
out_dir = fullfile(here, 'Model', 'Generated');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fprintf('\n=== ch5_gen_pendulum ===\n');
t0 = tic;

%% ------------------------------------------------- symbolic model, one source
% Note this calls ch5_pendulum -- the SAME function the simulation uses. There
% is no second, symbolic-only statement of the dynamics to fall out of sync.
x  = sym('x',  [8 1], 'real');
pv = sym('pv', [12 1], 'real');

d = ch5_pendulum(x, pv);
f = d.f;
g = d.g;

fprintf(' model built (%.1f s)\n', toc(t0));

%% ------------------------------------------------------------------ theta
% The controlled outputs are y = theta - theta_d, so y and theta share every
% derivative from the first on; theta_d never enters the generated code.
[th_lie, th_rd] = lie_stack(d.theta, f, g, x, 4);
fprintf(' theta: relative degree %d (expected 4)\n', th_rd);
assert(th_rd == 4, 'ch5_gen_pendulum:relDegTheta', ...
       'theta came out relative degree %d, not 4.', th_rd);

%% --------------------------------------------------------- end-effector height
[py_lie, py_rd] = lie_stack(d.py, f, g, x, 4);
fprintf(' py2:   relative degree %d (expected 4)\n', py_rd);
assert(py_rd == 4, 'ch5_gen_pendulum:relDegPy', ...
       'py2 came out relative degree %d, not 4.', py_rd);

%% ------------------------------------------------------------------- write
write_fn(fullfile(out_dir, 'ch5_pend_theta_lie.m'), ...
         {th_lie.d2, th_lie.d3, th_lie.drift, th_lie.in}, ...
         {'thd2', 'thd3', 'thd4_drift', 'thd4_in'}, x, pv, o.optimize, ...
         {'CH5_PEND_THETA_LIE  Exact time derivatives of theta, GENERATED CODE.', ...
          '', ...
          '  [thd2, thd3, thd4_drift, thd4_in] = ch5_pend_theta_lie(x, pv)', ...
          '', ...
          'For the 2-link pendulum with elastic actuators of Fig. 5.2:', ...
          '', ...
          '  thd2       = theta^(2)                    2x1', ...
          '  thd3       = theta^(3)                    2x1', ...
          '  thd4_drift = L_f^4 theta                  2x1', ...
          '  thd4_in    = L_g L_f^3 theta              2x2  (the decoupling matrix)', ...
          '', ...
          'so that  theta^(4) = thd4_drift + thd4_in * tau,  which is (5.21) with', ...
          'theta in place of B. thd4_in is the L_g L_f^(r-1) y of the (IO) row in', ...
          'every QP in this chapter; it equals D(theta)^-1 * k / Jm and is', ...
          'invertible wherever D is, i.e. everywhere.', ...
          '', ...
          'DO NOT EDIT. Regenerate with ch5_gen_pendulum.', ...
          'Checked against the live model in ch5_test_barrier.'});

write_fn(fullfile(out_dir, 'ch5_pend_py_lie.m'), ...
         {py_lie.d0, py_lie.d1, py_lie.d2, py_lie.d3, py_lie.drift, py_lie.in}, ...
         {'py', 'pyd1', 'pyd2', 'pyd3', 'pyd4_drift', 'pyd4_in'}, ...
         x, pv, o.optimize, ...
         {'CH5_PEND_PY_LIE  Exact time derivatives of the end-effector height, GENERATED.', ...
          '', ...
          '  [py, pyd1, pyd2, pyd3, pyd4_drift, pyd4_in] = ch5_pend_py_lie(x, pv)', ...
          '', ...
          'py2 of Fig. 5.2, world frame, up-positive, and its first four time', ...
          'derivatives along the closed-loop dynamics:', ...
          '', ...
          '  py .. pyd3   the eta_b stack of (5.10) once the constant p2min is', ...
          '               subtracted from py (a constant changes no derivative)', ...
          '  pyd4_drift = L_f^4 h                      1x1', ...
          '  pyd4_in    = L_g L_f^3 h                  1x2', ...
          '', ...
          'so that  h^(4) = pyd4_drift + pyd4_in * tau, the (5.21) split that makes', ...
          'the ECBF row affine in the control.', ...
          '', ...
          'DO NOT EDIT. Regenerate with ch5_gen_pendulum.', ...
          'Checked against the live model in ch5_test_barrier.'});

fprintf(' done in %.1f s -> %s\n\n', toc(t0), out_dir);

end

% ---------------------------------------------------------------------------
function [s, rd] = lie_stack(phi, f, g, x, r_expect)
%LIE_STACK  Build L_f^i phi for i = 0..r and the input coefficient at level r.
%
%   psi_0     = phi
%   psi_{i+1} = d(psi_i)/dx * f            <- pure drift, no control yet
%   L_g L_f^(r-1) phi = d(psi_{r-1})/dx * g
%
% The relative degree is the FIRST level at which the input coefficient is not
% identically zero, and it is measured rather than assumed: if the model is
% edited so that the control reaches the output sooner, this reports the new
% number and the caller's assert fires, instead of the generated code silently
% describing a system that no longer exists.

psi = phi;
s   = struct();
rd  = NaN;

for i = 1:r_expect
    Ji   = jacobian(psi, x);
    coef = Ji * g;

    if ~is_sym_zero(coef)
        rd = i;
        s.drift = Ji * f;
        s.in    = simplify(coef);      % small, and worth the tidy
        break;
    end

    psi = Ji * f;
    s.(sprintf('d%d', i)) = psi;
end

s.d0 = phi;

if isnan(rd)
    rd = r_expect + 1;      % control had still not appeared; caller asserts
end

end

% ---------------------------------------------------------------------------
function tf = is_sym_zero(e)
%IS_SYM_ZERO  Exact zero test, cheap path first.
%
% Below the relative degree the coefficient is STRUCTURALLY zero -- psi simply
% does not mention the states g touches -- so the trivial comparison settles it
% without calling simplify on a fourth-order trigonometric expression. simplify
% is only reached in the case where the answer is genuinely nonobvious, which
% for these two plants never happens.
tf = isequaln(e, sym(zeros(size(e))));
if ~tf
    tf = isequaln(simplify(e), sym(zeros(size(e))));
end
end

% ---------------------------------------------------------------------------
function write_fn(path, exprs, names, x, pv, do_opt, comments)
%WRITE_FN  matlabFunction wrapper with a real docstring and a size report.

fprintf(' writing %s ...', names{1});
t = tic;

matlabFunction(exprs{:}, ...
    'File',     path, ...
    'Vars',     {x, pv}, ...
    'Outputs',  names, ...
    'Optimize', do_opt, ...
    'Comments', comments);

info = dir(path);
fprintf(' %.1f kB (%.1f s)\n', info.bytes/1024, toc(t));

end
