function x_plus = rabbit_impact_map(x_minus,p)

q  = x_minus(1:7);
dq = x_minus(8:14);

D = D_matrix(q,p);
J = J_swing(q,p);   % <--- USE SWING FOOT JACOBIAN HERE

A = [D -J';
     J zeros(size(J,1))];

b = [D*dq;
     zeros(size(J,1),1)];

sol = A\b;

dq_plus = sol(1:7);

x_plus = [q; dq_plus];

end
