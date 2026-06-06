function [y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq)
% Maps the 7-DOF state to the 4 actuated joint outputs
% Rabbit q: [x, z, theta, q1, q2, q3, q4]
H_0 = [zeros(4, 3), eye(4)];
y = H_0 * q;
J_y = H_0;
Jdot_y_dq = zeros(4, 1);
end
