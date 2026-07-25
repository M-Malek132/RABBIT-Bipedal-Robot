function F = ch3_forces(t, X, alpha, p)
%CH3_FORCES  Recover torques, GRF and the Table 3.1 quantities along a trajectory.
%
%   F = ch3_forces(t, X, alpha, p)
%
% The ODE solver returns only states, so the torques and contact forces have
% to be recomputed afterwards.  Doing that HERE, through the same ch3_control
% entry point the simulation used, is what keeps the reported quantities
% honest: change p.controller and this function reports the new controller's
% forces automatically, with no second copy of the control law to drift out of
% sync.
%
% Inputs
%   t     : 1 x nt time grid
%   X     : 14 x nt state trajectory
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%
% Outputs
%   F : struct with columns aligned to t --
%         .u        nu x nt joint torques [Nm]
%         .lambda   2 x nt  ground reaction [Fx; Fz] world frame [N]
%         .mu_req   1 x nt  |Fx|/Fz, the friction coefficient DEMANDED
%         .y .ydot  ny x nt virtual constraint error and rate
%         .V        1 x nt  CLF value (NaN for non-CLF controllers)
%         .delta    1 x nt  CLF slack used (stage 8 only)
%         .s        1 x nt  phase
%         .sw_height 1 x nt swing foot height [m]
%         .int_u2   scalar  int ||u||^2 dt  (the cost numerator)
%         .torque_max, .Fz_min, .mu_max, .qp_infeasible summary scalars
%
% See also CH3_CONTROL, CH3_REPORT.

nt = numel(t);
nu = p.nu;
ny = p.ny;

F = struct();
F.u         = zeros(nu, nt);
F.lambda    = zeros(2,  nt);
F.mu_req    = zeros(1,  nt);
F.y         = zeros(ny, nt);
F.ydot      = zeros(ny, nt);
F.V         = nan(1,   nt);
F.delta     = zeros(1, nt);
F.s         = zeros(1, nt);
F.sw_height = zeros(1, nt);

n_infeas = 0;

for k = 1:nt
    x = X(:, k);
    [u, ci] = ch3_control(x, alpha, p);

    F.u(:, k)   = u;
    lam         = ci.aux.lam_drift + ci.aux.lam_in * u;
    F.lambda(:, k) = lam;

    % Friction coefficient DEMANDED by this contact force. Guarded: when Fz
    % passes through zero the ratio is unbounded and meaningless (the foot is
    % lifting off, not sliding), so it is reported as Inf rather than a large
    % finite number that would look like a merely-strict friction requirement.
    if lam(2) > 1e-6
        F.mu_req(k) = abs(lam(1)) / lam(2);
    else
        F.mu_req(k) = Inf;
    end

    F.y(:, k)    = ci.y;
    F.ydot(:, k) = ci.ydot;
    F.V(k)       = ci.V;
    F.delta(k)   = ci.delta;
    F.s(k)       = ci.o.s;
    if ~ci.qp_feasible, n_infeas = n_infeas + 1; end

    Psw = P_sw(x(1:p.nq));
    F.sw_height(k) = Psw(2);
end

% int ||u||^2 dt by the trapezoid rule on the solver's own grid
u2 = sum(F.u.^2, 1);
if nt > 1
    F.int_u2 = trapz(t, u2);
else
    F.int_u2 = 0;
end

F.torque_max     = max(abs(F.u(:)));
F.Fz_min         = min(F.lambda(2, :));
F.Fz_max         = max(F.lambda(2, :));
F.mu_max         = max(F.mu_req(isfinite(F.mu_req)));
if isempty(F.mu_max), F.mu_max = Inf; end
F.qp_infeasible  = n_infeas;
F.delta_max      = max(F.delta);

end
