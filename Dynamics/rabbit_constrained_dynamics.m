function [ddq,lambda] = rabbit_constrained_dynamics(q,dq,u)

% Dynamics terms
D = D(q);

C = C(q,dq);

G = G(q);

B = input_matrix();

% Contact Jacobian
J = J_st(q);

% Jdot*dq term
Jdotdq = Jdotdq_st(q,dq);

% Augmented system
A = [D   -J';
     J    zeros(size(J,1))];

b = [B*u - C - G;
     -Jdotdq];

% Solve system
x = A\b;

% Extract results
nq = length(q);

ddq = x(1:nq);

lambda = x(nq+1:end);

end
