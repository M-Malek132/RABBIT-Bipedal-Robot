function [ddq,lambda] = rabbit_constrained_dynamics(q,dq,u,foot_ref,alpha,beta)
%RABBIT_CONSTRAINED_DYNAMICS  Single-support forward dynamics (stance foot pinned).
%
%   [ddq,lambda] = rabbit_constrained_dynamics(q,dq,u)
%       Acceleration-level contact constraint only: J_st*ddq + Jdot*dq = 0.
%       This is a DOUBLE INTEGRATOR in the constraint, with no restoring term,
%       so any numerical residual integrates twice into POSITION drift -- the
%       stance foot slowly slides during a step (observed ~0.1-0.3 m).
%
%   [ddq,lambda] = rabbit_constrained_dynamics(q,dq,u,foot_ref,alpha,beta)
%       Adds BAUMGARTE STABILIZATION about the pinned location foot_ref (2x1,
%       the world-frame stance-foot position where the foot should stay, i.e.
%       P_st at the start of the step). The constraint becomes
%           c_ddot + 2*alpha*c_dot + beta^2*c = 0,   c = P_st(q) - foot_ref,
%       so position/velocity violations decay (critically damped when
%       alpha == beta) instead of accumulating. alpha,beta default to 20.
%
%   Passing foot_ref = [] (or omitting it) reproduces the original
%   acceleration-only behaviour, so existing callers are unaffected.

% Dynamics terms (renamed to avoid collision with function names)
M_mat = M(q);
C_vec = V([q;dq]);
G_vec = G(q);
B_mat = input_matrix();

% Contact Jacobian
J      = J_st(q);
Jdotdq = Jdotdq_st(q,dq);

% Baumgarte constraint-stabilization term (0 if no reference given)
if nargin >= 4 && ~isempty(foot_ref)
    if nargin < 5 || isempty(alpha), alpha = 20; end
    if nargin < 6 || isempty(beta),  beta  = 10; end
    c_pos = P_st(q) - foot_ref(:);   % position violation (foot vs pin)   2x1
    c_vel = J * dq;                  % velocity violation (foot velocity) 2x1
    stab  = 2*alpha*c_vel + beta^2*c_pos;
else
    stab  = zeros(size(J,1),1);
end

% Augmented (KKT) system.  Bottom row enforces:
%   J*ddq = -Jdot*dq - (2*alpha*c_dot + beta^2*c)
A = [M_mat   -J';
     J        zeros(size(J,1))];

b = [B_mat*u - C_vec - G_vec;
     -Jdotdq - stab];

% Solve
x = A\b;

nq     = length(q);
ddq    = x(1:nq);
lambda = x(nq+1:end);

end
