function u = rabbit_fl_controller(t, x, p)
    % Extract joint positions and velocities
    q  = x(1:7);
    dq = x(8:14);

    % Get constrained dynamic matrices: 
    % Dc*ddq + Hc = Bc*u 
    [Dc, Hc, Bc] = rabbit_constrained_dynamics(q, dq, p);

    % Get outputs, their derivatives, and Jacobians
    % y = h(q), dy = Hq*dq
    [y, dy, Hq, dHq] = rabbit_outputs(q, dq, p);

    % Decoupling matrix A and drift vector b
    % \ddot{y} = A*u + b
    A = Hq * (Dc \ Bc);
    b = -Hq * (Dc \ Hc) + dHq * dq;

    % PD gains for the virtual control input v
    % Kp=50 and Kd=15 provide a slightly overdamped response
    % (critical damping is Kd = 2*sqrt(Kp) ~ 14.1)
    Kp = 50 * eye(4);
    Kd = 15 * eye(4);

    % Virtual control law
    v = -Kp * y - Kd * dy;

    % Actual control input (Feedback Linearization)
    % Solves Au + b = v for u
    u = A \ (v - b);
end
