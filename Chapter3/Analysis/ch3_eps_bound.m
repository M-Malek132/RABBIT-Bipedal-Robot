function B = ch3_eps_bound(file)
%CH3_EPS_BOUND  The largest eps the RES-CLF can certify through an impact.
%
%   B = ch3_eps_bound()        on Results/ch3_gait_posture_195.mat
%   B = ch3_eps_bound(file)    on any gait file (z, p)
%
% Chapter 3 CHOOSES eps (0.5, then 0.20 in Chapter 4) and explains afterwards
% why the min-norm CLF-QP falls at 0.5. This DERIVES an upper bound on eps
% from the same theorem, so the choice has a number to be compared against.
%
% THE ARGUMENT. Between impacts the RES-CLF condition gives
%       V(T-) <= exp(-c3 T / eps) V(0+),
% and the impact maps the transverse state eta = [y; ydot] linearly near the
% fixed point, eta+ = J eta-, so
%       V(eta+) <= gamma2(eps) V(eta-),
%       gamma2   = max_eta V(J eta) / V(eta) = lambda_max(J' P_eps J, P_eps).
% One whole step therefore multiplies V by at most
%       rho_V(eps) = gamma2(eps) exp(-c3 T / eps),
% and the RES-CLF certifies the hybrid loop only where rho_V < 1. The
% conservative form ||J||^2 cond(P_eps) exp(-c3 T/eps) is reported too: it is
% the bound one writes down without solving the generalized eigenproblem.
% Both hold in CONTINUOUS time and at the fixed point; sampling and distance
% from the orbit only make things worse, so eps_max is an upper bound on what
% can work, and the gap to what works in simulation is itself a measurement.
%
% HOW J IS MEASURED. At the pre-impact state of the gait's own fixed point,
% states are built with the zero-dynamics coordinates held (theta, thetadot
% at their fixed-point values), the stance foot where it was, and a chosen
% eta: q_act = yd(s) + y, the torso from c_theta q = theta, the floating base
% from the foot, and qdot from [Jy; c_theta; J_st] qdot = [ydot; thetadot; 0].
% Each goes through ch3_impact, and J is the central difference of the eta
% that comes out. Only the reset is linearized (the guard-timing shift of a
% full Poincare linearization is not); it is the review's suggested bound.
%
% Output: struct B, and Results/reruns/ch3_eps_bound/eps_bound.log (seconds).
%
% See also CH3_RES_CLF, CH3_IMPACT, CH3_ZD_POINT, CH3_POINCARE.

ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1 || isempty(file)
    file = fullfile(ROOT, 'Results', 'ch3_gait_posture_195.mat');
end
OUT = fullfile(ROOT, 'Results', 'reruns', 'ch3_eps_bound');
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'eps_bound.log');
if exist(logf, 'file'), delete(logf); end

S = load(file);
p = ch3_upgrade_params(S.p);
E = ch3_col_eval(S.z, p);
p = E.p;
[X, T, alpha] = ch3_col_unpack(S.z, p);
xm = X(:, end);                                  % pre-impact, on Z, at the guard

[J, eta_m, eta_p] = impact_eta_jacobian(xm, alpha, p);
logln(logf, '=== ch3_eps_bound | %s | T %.4f s | %s', file, T, datestr(now)); %#ok<TNOW1,DATST>
logln(logf, ['impact on eta at the fixed point: ||J||_2 = %.3f, spectral radius %.3f ' ...
             '(|eta-| %.1e, |eta+| %.1e on the orbit)'], norm(J), max(abs(eig(J))), ...
      norm(eta_m), norm(eta_p));

B = struct('file', file, 'T', T, 'J', J, 'normJ', norm(J));
eps_grid = logspace(log10(0.01), log10(1.0), 600);

for cons = {'care', 'lyap'}
    % P and c3 do not depend on eps; P_eps = I_eps P I_eps, I_eps = diag(I/eps, I),
    % exactly as ch3_res_clf builds it (done here to keep its cache small).
    pc = p; pc.clf_construction = cons{1}; pc.eps = 1;
    clf1 = ch3_res_clf(pc);
    ny = p.ny;
    g_t = zeros(size(eps_grid)); g_c = g_t; dec = g_t;
    for i = 1:numel(eps_grid)
        Ie = blkdiag(eye(ny) / eps_grid(i), eye(ny));
        Pe = Ie * clf1.P * Ie;
        Pe = (Pe + Pe.') / 2;
        g_t(i) = max(real(eig(J.' * Pe * J, Pe)));      % generalized, tight
        g_c(i) = norm(J)^2 * cond(Pe);                  % conservative
        dec(i) = exp(-clf1.c3 * T / eps_grid(i));
    end
    rho_t = g_t .* dec;
    rho_c = g_c .* dec;
    et = eps_max(eps_grid, rho_t);
    ec = eps_max(eps_grid, rho_c);
    at = @(e) interp1(eps_grid, rho_t, e);
    ac = @(e) interp1(eps_grid, rho_c, e);
    logln(logf, ['%-4s  c3 = %.4f | eps_max: %.4f (generalized eigenvalue), %.4f ' ...
                 '(||J||^2 cond P_eps)'], upper(cons{1}), clf1.c3, et, ec);
    logln(logf, ['      per-step V multiplier at eps 0.50: %.3g (tight), %.3g (conservative); ' ...
                 'at 0.20: %.3g, %.3g; at 0.10: %.3g, %.3g'], ...
          at(0.5), ac(0.5), at(0.2), ac(0.2), at(0.1), ac(0.1));
    B.(cons{1}) = struct('eps', eps_grid, 'gamma2', g_t, 'gamma2_cons', g_c, ...
                         'decay', dec, 'rho', rho_t, 'rho_cons', rho_c, ...
                         'eps_max', et, 'eps_max_cons', ec, 'c3', clf1.c3);
end
save(fullfile(OUT, 'eps_bound.mat'), 'B');
logln(logf, '=== EPS_BOUND_DONE');
end

% ---------------------------------------------------------------------------
function e = eps_max(grid, rho)
%EPS_MAX  Largest eps below which rho < 1 on the whole grid (NaN: never).
bad = find(rho >= 1, 1, 'first');
if isempty(bad),  e = grid(end); return; end
if bad == 1,      e = NaN;       return; end
% linear interpolation of the crossing in log(rho)
i0 = bad - 1;
t  = (0 - log(rho(i0))) / (log(rho(bad)) - log(rho(i0)));
e  = grid(i0) + t * (grid(bad) - grid(i0));
end

% ---------------------------------------------------------------------------
function [J, eta_m, eta_p] = impact_eta_jacobian(xm, alpha, p)
nq = p.nq; ny = p.ny;
c  = p.c_theta(:).';
q0 = xm(1:nq); dq0 = xm(nq+1:end);
th0 = c * q0; thd0 = c * dq0;
[y0, yd0] = ch3_outputs(xm, alpha, p);
eta_m = [y0; yd0];
xp0 = ch3_impact(xm, p);
[yp, ydp] = ch3_outputs(xp0, alpha, p);
eta_p = [yp; ydp];

J = zeros(2*ny);
for i = 1:2*ny
    h  = 1e-6;
    ep = zeros(2*ny, 1); ep(i) = h;
    x1 = state_from_eta(q0, th0, thd0, eta_m + ep, alpha, p);
    x2 = state_from_eta(q0, th0, thd0, eta_m - ep, alpha, p);
    [y1, yd1] = ch3_outputs(ch3_impact(x1, p), alpha, p);
    [y2, yd2] = ch3_outputs(ch3_impact(x2, p), alpha, p);
    J(:, i) = ([y1; yd1] - [y2; yd2]) / (2*h);
end
end

% ---------------------------------------------------------------------------
function x = state_from_eta(q0, th, thd, eta_t, alpha, p)
%STATE_FROM_ETA  The state with phase (th, thd), transverse state eta_t, and
% the stance foot where it is at q0 (at rest on the ground).
nq = p.nq; ny = p.ny;
c   = p.c_theta(:).';
dth = p.theta_plus - p.theta_minus;
s   = (th - p.theta_minus) / dth;
[yd, dyd] = ch3_yd(alpha, s, p);

q = q0;
q(p.iact) = yd + eta_t(1:ny);                   % y = q_act - yd(s)
i_free = find(c ~= 0 & ~ismember(1:nq, p.iact(:).'));
q(i_free) = (th - c(p.iact) * q(p.iact)) / c(i_free);
Jb = J_st(q);
q(1:2) = q(1:2) - Jb(:, 1:2) \ (P_st(q) - P_st(q0));   % affine: exact

Jy = p.H - dyd * (c / dth);
A  = [Jy; c; J_st(q)];
qd = A \ [eta_t(ny+1:end); thd; zeros(2, 1)];
x  = [q; qd];
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
