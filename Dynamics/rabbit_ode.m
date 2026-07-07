function dx = rabbit_ode(t, x, controller)


% split state
q  = x(1:7);
dq = x(8:14);

% control input
if isempty(controller)
    u = zeros(4,1);           % e.g., 4 joint torques
else
    u = controller(x);
end

% call constrained dynamics
ddq = rabbit_constrained_dynamics(q, dq, u);

% assemble state derivative
dx = [dq; ddq];
end
