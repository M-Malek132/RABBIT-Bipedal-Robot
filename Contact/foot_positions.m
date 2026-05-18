function [p_stance,p_swing] = foot_positions(q,p_robot)

L1 = p_robot.l1;
L2 = p_robot.l2;

x  = q(1);
z  = q(2);

qt = q(3);

q1 = q(4);
q2 = q(5);

q3 = q(6);
q4 = q(7);

% Absolute link angles
th1 = qt + q1;
th2 = th1 + q2;

th3 = qt + q3;
th4 = th3 + q4;

% Stance foot
p_stance = [
    x - L1*sin(th1) - L2*sin(th2);
    z - L1*cos(th1) - L2*cos(th2)
];

% Swing foot
p_swing = [
    x - L1*sin(th3) - L2*sin(th4);
    z - L1*cos(th3) - L2*cos(th4)
];

end
