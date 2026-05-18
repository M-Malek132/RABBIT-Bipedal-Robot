function [ddq,lambda] = rabbit_constrained_dynamics(q,dq,u,param)

% Dynamics terms
D = D_matrix(q,param);

C = C_vector(q,dq,param);

G = G_vector(q,param);

B = input_matrix();

% Contact Jacobian
J = J_stance(q,param);

% Jdot*dq term
Jdotdq = Jdotdq_stance(q,dq,param);

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
