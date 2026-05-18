function dx = rabbit_ode(t, x, param, controller)

% convert struct → numeric vector for use in low-level functions
p = packParameters(param);

% split state
q  = x(1:7);
dq = x(8:14);

% control input
if isempty(controller)
    u = zeros(4,1);           % e.g., 4 joint torques
else
    u = controller(t, x, param);
end

% call constrained dynamics
ddq = rabbit_constrained_dynamics(q, dq, u, p);

% assemble state derivative
dx = [dq; ddq];
end
