function [ddq,lambda] = rabbit_constrained_dynamics(q,dq,u)

% Dynamics terms
% Rename variables to avoid collision with function names
M_mat = M(q);      % Was D = M(q)
C_vec = V([q;dq]); % Was C = V([q;dq])
G_vec = G(q);      % Was G = G(q)
B_mat = input_matrix();

% Contact Jacobian
J = J_st(q);
Jdotdq = Jdotdq_st(q,dq);

% Augmented system
% Use the new variable names
A = [M_mat   -J';
     J        zeros(size(J,1))];

b = [B_mat*u - C_vec - G_vec;
     -Jdotdq];

% Solve system
x = A\b;

% Extract results
nq = length(q);
ddq = x(1:nq);
lambda = x(nq+1:end);

end
