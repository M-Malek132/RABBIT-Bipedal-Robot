function [p_stance,p_swing,p_hip,p_knee_s,p_knee_sw,p_torso] = ...
    rabbit_kinematics(q,param)

px = q(1);
pz = q(2);

qt = q(3);

q1 = q(4);
q2 = q(5);

q3 = q(6);
q4 = q(7);

L1 = param(4);
L2 = param(5);
Lt = param(6);

%% Hip

p_hip = [px; pz];

%% Stance leg

p_knee_s = p_hip - [
    L1*sin(qt + q1);
    L1*cos(qt + q1)
];

p_stance = p_knee_s - [
    L2*sin(qt + q1 + q2);
    L2*cos(qt + q1 + q2)
];

%% Swing leg

p_knee_sw = p_hip - [
    L1*sin(qt + q3);
    L1*cos(qt + q3)
];

p_swing = p_knee_sw - [
    L2*sin(qt + q3 + q4);
    L2*cos(qt + q3 + q4)
];

%% Torso

p_torso = p_hip + [
    Lt*sin(qt);
    Lt*cos(qt)
];

end
