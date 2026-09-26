function D = ch3_stride_licq(src, cut)
%CH3_STRIDE_LICQ  Why the stride will not come in, tested at the stuck iterate.
%
%   D = ch3_stride_licq()             the iterate the stride march stalled at
%   D = ch3_stride_licq(file)         any gait or iterate .mat (z | z_try | z_opt, p)
%   D = ch3_stride_licq(file, cut)    first-order stride reduction asked [m] (1e-3)
%
% The stride march (ch3_stride_march) could not shorten the stride by even a
% millimetre: the row-19 violation stayed constant, the equalities held and
% SQP's step collapsed to 1e-8. Four explanations were proposed and falsified
% (see docs/HANDOFF.md). This tests two more, AT the point where it stalled,
% before anyone re-solves anything.
%
% 1. GEOMETRY. theta = qt + q1 + q2/2 is the hip-to-foot direction (0.5 m
%    links), both feet are on the ground at the strike, and the landing leg is
%    at theta_minus after relabelling, so L_step = h_imp (tan theta_plus -
%    tan theta_minus), h_imp the hip height at impact. Printed both ways, with
%    the floor the hip band (row 8) puts under the stride.
%
% 2. LICQ (the reviewers' hypothesis). The Jacobian of the ACTIVE constraints
%    -- all equalities, the active inequalities, the active bounds -- and its
%    smallest singular values with and without the stride row, plus the part
%    of grad L orthogonal to the others' row space. The aggregated rows
%    (max/min over nodes: torque, friction, Fz, hip band, clearance, theta-dot)
%    are DISAGGREGATED here, one row per node and midpoint (and per joint and
%    side), so a kink of the max() is not mistaken for rank loss.
%
% 3. THE LINEARIZED PROBLEM. Is there ANY first-order direction d that keeps
%    every equality (J_eq d = 0) and every active inequality (J_act d <= 0),
%    respects the bounds, and shortens the stride by `cut`? Solved as an
%    elastic LP (minimize the shortfall s). s = 0: such a direction exists and
%    the stall is not the constraint geometry. s > 0: none does, and the LP's
%    multipliers name the constraints that block it -- the answer to "what is
%    the obstacle", row by row.
%
% 4. The same three with theta_minus, theta_plus free (p.free_theta), i.e. with
%    the two columns the free-theta transcription adds -- what freeing them
%    would change at this very point.
%
% Cost: a central-difference Jacobian, two per case: ~3500 evaluations of
% ch3_col_constraints at N = 61, about 5-10 minutes. Output:
% Results/reruns/stride_licq/ (licq.log and licq.mat).
%
% See also CH3_STRIDE_MARCH, CH3_STRIDE_FREE_THETA, CH3_COL_CONSTRAINTS.

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
if nargin < 1 || isempty(src)
    src = fullfile(ROOT, 'Results', 'reruns', 'stride_march', 'stride_L400.mat');
end
if nargin < 2 || isempty(cut), cut = 1e-3; end
OUT = fullfile(ROOT, 'Results', 'reruns', 'stride_licq');
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'licq.log');
if exist(logf, 'file'), delete(logf); end        % diary-style logs append

S = load(src);
z = [];
for c = {'z_try', 'z', 'z_opt'}
    if isfield(S, c{1}), z = S.(c{1})(:); break; end
end
if isempty(z) || ~isfield(S, 'p')
    error('ch3_stride_licq:contents', '"%s" needs z (or z_try, z_opt) and p.', src);
end
p = ch3_upgrade_params(S.p);
p.free_theta = false;
z = z(1:end - 2*(mod(numel(z) - 1 - p.ny*p.n_ctrl, p.nx) ~= 0));   % drop a theta tail

logln(logf, '=== ch3_stride_licq | %s | cut %.1e m | %s', src, cut, datestr(now)); %#ok<TNOW1,DATST>

%% 1. geometry ---------------------------------------------------------------
E  = ch3_col_eval(z, p);
nq = p.nq;
c  = p.c_theta(:).';
th1 = c * E.X(1:nq, 1);
thN = c * E.X(1:nq, E.N);
swN = P_sw(E.X(1:nq, E.N));
h_imp = -E.X(2, E.N) - swN(2);                 % hip above the landing foot
L_pred = h_imp * (tan(thN) - tan(th1));
logln(logf, ['GEOMETRY  L_step %.6f m | h_imp (tan th+ - tan th-) = %.6f x %.6f = %.6f m ' ...
             '| difference %.1e'], E.L_step, h_imp, tan(thN) - tan(th1), L_pred, E.L_step - L_pred);
if p.limits.enable.height
    h_lo = p.limits.hip_h - p.limits.hip_h_tol;
    logln(logf, ['          hip band on: h_imp >= %.3f m, so with theta fixed the stride ' ...
                 'cannot go below %.4f m'], h_lo, h_lo * (tan(thN) - tan(th1)));
end
D = struct('src', src, 'cut', cut, 'L', E.L_step, 'h_imp', h_imp, ...
           'L_pred', L_pred, 'theta', [th1 thN]);

%% 2-4. fixed theta, then free theta --------------------------------------
cases = {'fixed theta', p, z};
pf = p; pf.free_theta = true;
cases(2, :) = {'free theta', pf, ch3_col_theta_augment(z, pf)};

for k = 1:2
    [label, pk, zk] = cases{k, :};
    logln(logf, ' ');
    logln(logf, '--- %s: %d variables ---', label, numel(zk));
    t0 = tic;
    R = analyse(zk, pk, cut, logf);
    R.seconds = toc(t0);
    logln(logf, '    (%.0f s)', R.seconds);
    D.(matlab.lang.makeValidName(label)) = R;
end

save(fullfile(OUT, 'licq.mat'), 'D');
logln(logf, '=== STRIDE_LICQ_DONE');
fprintf('STRIDE_LICQ_DONE\n');
end

% ===========================================================================
function R = analyse(z, p, cut, logf)
[F0, names, kind] = rows_at(z, p);
n  = numel(z);
E0 = ch3_col_eval(z, p);
[lb, ub] = ch3_col_bounds(p, E0.N);

% --- central-difference Jacobian of every row -----------------------------
J = zeros(numel(F0), n);
for j = 1:n
    h  = 1e-6 * max(1, abs(z(j)));
    zp = z; zp(j) = zp(j) + h;
    zm = z; zm(j) = zm(j) - h;
    J(:, j) = (rows_at(zp, p) - rows_at(zm, p)) / (2*h);
end

is_eq  = strcmp(kind, 'eq');
is_L   = strcmp(kind, 'L');
is_in  = strcmp(kind, 'in');
gradL  = J(is_L, :);

% equalities that are constant (a gate held at 0, e.g. NEC1 off) carry no
% information and would only add zero rows
Jeq = J(is_eq, :);
keep_eq = any(abs(Jeq) > 0, 2);
Jeq = Jeq(keep_eq, :);

% active inequalities: within 1e-6 of their bound, in the units of the row
tol_act = 1e-6 * max(1, abs(F0));
act = is_in & (F0 >= -tol_act);
near = is_in & (F0 >= -1e-3 * max(1, abs(F0)));
Jact = J(act, :);
act_names = names(act);

% active bounds, as rows of the identity
at_lo = abs(z - lb) <= 1e-9 * max(1, abs(lb));
at_hi = abs(ub - z) <= 1e-9 * max(1, abs(ub));
Ib = eye(n);
Jb = [-Ib(at_lo, :); Ib(at_hi, :)];
b_names = [arrayfun(@(j) sprintf('bound: z(%d) at lower', j), find(at_lo), 'UniformOutput', false); ...
           arrayfun(@(j) sprintf('bound: z(%d) at upper', j), find(at_hi), 'UniformOutput', false)];

logln(logf, ['    active set: %d equalities (of %d), %d inequality rows active ' ...
             '(%d within 1e-3), %d bounds'], size(Jeq, 1), nnz(is_eq), nnz(act), ...
      nnz(near), size(Jb, 1));
summarize_types(logf, act_names);

% --- 2. LICQ ---------------------------------------------------------------
A  = [Jeq; Jact; Jb];
sA = svd(A);
sB = svd([A; gradL]);
U  = orth(A.');                                % basis of the row space of A
r_perp = norm(gradL.' - U * (U.' * gradL.')) / norm(gradL);
R = struct('n_rows', size(A, 1), 'n_vars', n, ...
           'sig_min', sA(end) / sA(1), 'sig_min_with_L', sB(end) / sB(1), ...
           'rank', rank(A), 'rank_with_L', rank([A; gradL]), 'r_perp', r_perp);
logln(logf, ['    LICQ: %d active rows on %d variables; rank %d -> %d with the stride row; ' ...
             'smallest singular value (relative) %.2e -> %.2e'], size(A, 1), n, ...
      R.rank, R.rank_with_L, R.sig_min, R.sig_min_with_L);
logln(logf, ['          |grad L off the active row space| / |grad L| = %.2e  ' ...
             '(0 = the stride row is a combination of the active ones)'], r_perp);

% --- 3. the linearized problem: can the stride come in by `cut`? ---------
% variables [d; s], minimize s
%   Jeq d = 0,  Jact d <= 0,  gradL d - s <= -cut,  lb-z <= d <= ub-z,
%   |d| <= 1 (a trust box, generous at these scales),  s >= 0
nv   = n + 1;
f    = [zeros(n, 1); 1];
Aeq  = [Jeq, zeros(size(Jeq, 1), 1)];
beq  = zeros(size(Jeq, 1), 1);
Ain  = [Jact, zeros(size(Jact, 1), 1); gradL, -1];
bin  = [zeros(size(Jact, 1), 1); -cut];
lbd  = [max(lb - z, -1); 0];
ubd  = [min(ub - z,  1); Inf];
opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');
[x, fval, flag, ~, lam] = linprog(f, Ain, bin, Aeq, beq, lbd, ubd, opts);
R.lp_flag = flag;
R.shortfall = fval;
if flag ~= 1
    logln(logf, '    LINEARIZED: linprog exit %d, no answer', flag);
    return;
end
if fval <= 1e-9 * cut
    dL = gradL * x(1:n);
    logln(logf, ['    LINEARIZED: FEASIBLE -- a first-order direction shortens the stride ' ...
                 'by %.2e m keeping every active constraint (|d|_inf %.2e).'], -dL, max(abs(x(1:n))));
    logln(logf, ['          So the constraint geometry does not forbid it here; the stall ' ...
                 'is curvature, the non-smooth max() rows or the line search.']);
    % which inequalities does that direction move along the boundary of?
    tight = find(abs(lam.ineqlin(1:end-1)) > 0);
    R.tight = act_names(tight);
else
    logln(logf, ['    LINEARIZED: INFEASIBLE -- no first-order direction shortens the stride ' ...
                 'by %.1e m; best shortfall %.3e m.'], cut, fval);
    % blocking rows: positive multipliers on the active inequalities and bounds
    mu_in = lam.ineqlin(1:end-1);
    [mu_s, order] = sort(mu_in, 'descend');
    top = order(mu_s > 1e-9 * max(1, max(mu_s)));
    logln(logf, '          blocking inequality rows (multiplier, name):');
    for i = 1:min(15, numel(top))
        logln(logf, '            %10.3e  %s', mu_in(top(i)), act_names{top(i)});
    end
    mb = [lam.lower(1:n); lam.upper(1:n)];
    bn = [arrayfun(@(j) sprintf('box: z(%d) lower', j), (1:n)', 'UniformOutput', false); ...
          arrayfun(@(j) sprintf('box: z(%d) upper', j), (1:n)', 'UniformOutput', false)];
    [mbs, ob] = sort(mb, 'descend');
    for i = 1:min(5, nnz(mbs > 1e-9))
        logln(logf, '            %10.3e  %s', mbs(i), bn{ob(i)});
    end
    R.blocking = act_names(top(1:min(15, numel(top))));
    R.blocking_mu = mu_in(top(1:min(15, numel(top))));
end
R.bound_names = b_names;
end

% ===========================================================================
function [F, names, kind] = rows_at(z, p)
%ROWS_AT  Every constraint as a smooth row: equalities, disaggregated
% inequalities (value <= 0 is feasible), and the stride. The names and kinds
% are only built when asked for.
[c, ceq] = ch3_col_constraints(z, p);
E  = ch3_col_eval(z, p);
p  = E.p;
N  = E.N;
L  = p.limits;
en = L.enable;
want = nargout > 1;

vals = {}; nm = {};

if en.torque
    U = [E.u, E.um];
    for j = 1:p.nu
        add( U(j,:) - L.u_max, per_node(sprintf('torque u%d upper', j), N, N-1));
        add(-U(j,:) - L.u_max, per_node(sprintf('torque u%d lower', j), N, N-1));
    end
end
lam = [E.lam, E.lamm];
if en.friction
    add( lam(1,:) - L.mu_s * lam(2,:), per_node('friction +Fx', N, N-1));
    add(-lam(1,:) - L.mu_s * lam(2,:), per_node('friction -Fx', N, N-1));
end
if en.grf
    add(L.Fz_min - lam(2,:), per_node('Fz floor', N, N-1));
end
if en.height
    hh = -[E.X(2,:), E.xm(2,:)];
    add( (hh - L.hip_h) - L.hip_h_tol, per_node('hip band top', N, N-1));
    add(-(hh - L.hip_h) - L.hip_h_tol, per_node('hip band bottom', N, N-1));
end
if en.swing_clear && N > 3
    hs = [E.sw_h(2:N-1), E.sw_hm(2:N-2)];
    add(L.sw_clear_min - hs, @(v) [arrayfun(@(i) sprintf('swing clear node %d', i), (2:N-1)', 'UniformOutput', false); ...
                                    arrayfun(@(i) sprintf('swing clear mid %d', i), (2:N-2)', 'UniformOutput', false)]);
end
if en.phase_mono
    add(L.thetadot_min - [E.thd, E.thdm], per_node('theta-dot floor', N, N-1));
end
if isfield(en, 'clearance_max') && en.clearance_max
    add([E.sw_h, E.sw_hm] - L.clearance_max, per_node('swing ceiling', N, N-1));
end
% the scalar rows, as ch3_col_constraints computes them (gated-off rows are
% held at -1 and can never be active)
scal = [1 3 7 10 11 12 13 14 15 17];
lbl  = {'row 1 mid-step clearance', 'row 3 stride floor', 'row 7 impulse', ...
        'row 10 strike rate', 'row 11 lift-off', 'row 12 Iz >= 0', ...
        'row 13 impact cone', 'row 14 NEC4 existence', 'row 15 NEC5 stability', ...
        'row 17 decoupling'};
for i = 1:numel(scal)
    add(c(scal(i)), @(v) lbl(i));
end
ineq = vertcat(vals{:});

F = [ceq(:); ineq; E.L_step];
if want
    in_names = vertcat(nm{:});
    names = [arrayfun(@(i) sprintf('equality %d', i), (1:numel(ceq))', 'UniformOutput', false); ...
             in_names; {'stride L'}];
    kind = [repmat({'eq'}, numel(ceq), 1); repmat({'in'}, numel(ineq), 1); {'L'}];
end

    % ---- nested helpers (share vals, nm and want with rows_at) -----------
    function add(v, label)
        vals{end+1} = v(:);
        if want, nm{end+1} = label(v); end
    end
    function lab = per_node(prefix, n1, n2)
        % names for [nodes 1..n1, midpoints 1..n2]
        lab = @(v) [arrayfun(@(i) sprintf('%s node %d', prefix, i), (1:n1)', 'UniformOutput', false); ...
                    arrayfun(@(i) sprintf('%s mid %d', prefix, i), (1:n2)', 'UniformOutput', false)];
    end
end

% ---------------------------------------------------------------------------
function summarize_types(logf, act_names)
if isempty(act_names), logln(logf, '    (no inequality active)'); return; end
stem = regexprep(act_names, ' (node|mid) \d+$', '');
[u, ~, g] = unique(stem);
cnt = accumarray(g, 1);
[cnt, o] = sort(cnt, 'descend');
logln(logf, '    active inequalities by type:');
for i = 1:numel(u)
    logln(logf, '      %4d  %s', cnt(i), u{o(i)});
end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
