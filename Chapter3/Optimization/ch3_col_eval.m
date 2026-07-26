function E = ch3_col_eval(z, p)
%CH3_COL_EVAL  Evaluate everything the cost and the constraints need, once.
%
%   E = ch3_col_eval(z, p)
%
% fmincon calls the objective and the nonlinear constraints back to back with
% the SAME z -- and does so twice per variable during finite differencing.
% Evaluating the nodes in both places would double the cost of the entire
% solve, so the work happens here and both callers read the cached result.
% The cache holds the last z only, which is all the access pattern needs.
%
% HERMITE-SIMPSON.  For each interval, with h = T/(N-1):
%
%       x_mid  = (x_k + x_k+1)/2 + (h/8)(f_k - f_k+1)
%       defect = x_k+1 - x_k - (h/6)(f_k + 4 f_mid + f_k+1)
%
% Driving the defects to zero IS the dynamics constraint -- there is no ODE
% solver anywhere in the optimizer loop.  Two consequences that matter:
% the stance foot is pinned at every node rather than only integrated, so the
% contact cannot drift during a step; and the same Simpson weights that define
% the defects also integrate the cost, so cost and dynamics use one consistent
% quadrature.
%
% Output E has fields:
%   X T alpha N h        unpacked decision variables
%   f  (14 x N)          node derivatives      u  (4 x N)      node torques
%   fm (14 x N-1)        midpoint derivatives  um (4 x N-1)    midpoint torques
%   xm (14 x N-1)        midpoint states
%   defect (14 x N-1)
%   lam (2 x N), lamm (2 x N-1)   ground reaction at nodes / midpoints
%   int_u2               int ||u||^2 dt by Simpson
%   L_step               step length [m]
%   sw_h (1 x N)         swing-foot height at each node
%   x_next, impulse      result of applying Delta at the final node
%
% See also CH3_COL_COST, CH3_COL_CONSTRAINTS, CH3_COL_DYNAMICS.

persistent key val
if ~isempty(key) && numel(key) == numel(z) && isequal(key, z)
    E = val;
    return;
end

[X, T, alpha] = ch3_col_unpack(z, p);
N  = size(X, 2);
h  = T / (N - 1);
nq = p.nq;

f   = zeros(p.nx, N);
u   = zeros(p.nu, N);
lam = zeros(2,    N);
sw_h = zeros(1,   N);

for k = 1:N
    [f(:,k), u(:,k), aux] = ch3_col_dynamics(X(:,k), alpha, p);
    lam(:,k) = aux.lam_drift + aux.lam_in * u(:,k);
    Psw = P_sw(X(1:nq, k));
    sw_h(k) = Psw(2);
end

xm     = zeros(p.nx, N-1);
fm     = zeros(p.nx, N-1);
um     = zeros(p.nu, N-1);
lamm   = zeros(2,    N-1);
defect = zeros(p.nx, N-1);

for k = 1:N-1
    xm(:,k) = 0.5*(X(:,k) + X(:,k+1)) + (h/8)*(f(:,k) - f(:,k+1));
    [fm(:,k), um(:,k), auxm] = ch3_col_dynamics(xm(:,k), alpha, p);
    lamm(:,k) = auxm.lam_drift + auxm.lam_in * um(:,k);

    defect(:,k) = X(:,k+1) - X(:,k) - (h/6)*(f(:,k) + 4*fm(:,k) + f(:,k+1));
end

% Simpson quadrature of ||u||^2, same weights as the defects
u2  = sum(u.^2,  1);
um2 = sum(um.^2, 1);
int_u2 = 0;
for k = 1:N-1
    int_u2 = int_u2 + (h/6)*(u2(k) + 4*um2(k) + u2(k+1));
end

% step length: stance foot at the start -> swing foot at the strike
foot_st_0  = P_st(X(1:nq, 1));
foot_sw_N  = P_sw(X(1:nq, N));
L_step     = foot_sw_N(1) - foot_st_0(1);

[x_next, impulse] = ch3_impact(X(:,N), p);

E = struct('X', X, 'T', T, 'alpha', alpha, 'N', N, 'h', h, ...
           'f', f, 'u', u, 'fm', fm, 'um', um, 'xm', xm, ...
           'defect', defect, 'lam', lam, 'lamm', lamm, ...
           'int_u2', int_u2, 'L_step', L_step, 'sw_h', sw_h, ...
           'x_next', x_next, 'impulse', impulse);

key = z;
val = E;

end
