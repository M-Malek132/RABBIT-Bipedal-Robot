clc
clear

%% parameters
p = parameters();
par = packParameters(p);

%% random configuration
q  = randn(7,1);
dq = randn(7,1);

eps = 1e-6;

%% stance foot function (same formula you used)
stanceFoot = @(q,par)[
    q(1) - p.l1*sin(q(3)+q(4)) - p.l2*sin(q(3)+q(4)+q(5));
    q(2) - p.l1*cos(q(3)+q(4)) - p.l2*cos(q(3)+q(4)+q(5))
];

%% compute position
p0 = stanceFoot(q,par);

%% finite difference velocity
p1 = stanceFoot(q + eps*dq,par);
fd = (p1 - p0)/eps;

%% Jacobian prediction
J = J_stance(q,par);
jac = J*dq;

%% display
disp('Finite difference velocity:')
disp(fd)

disp('Jacobian velocity:')
disp(jac)

disp('Error:')
disp(fd - jac)
