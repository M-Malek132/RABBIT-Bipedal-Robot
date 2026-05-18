clc
clear
close all

%% =========================
%% Generalized coordinates
%% =========================

syms q  [7 1] real
syms dq [7 1] real

%% =========================
%% Symbolic parameters
%% =========================

syms mT m1 m2 l1 l2 lt g real
syms I1 I2 IT real

param = [m1 m2 mT l1 l2 lt I1 I2 IT g];

%% =========================
%% Coordinates (stance leg = q4, q5)
%% =========================

px = q(1);
pz = q(2);

% stance leg
q_st_hip  = q(4);
q_st_knee = q(5);

% swing leg
q_sw_hip  = q(6);
q_sw_knee = q(7);

%% =========================
%% Stance foot position
%% =========================

p_st = [
    px - l1*sin(q(3) +q_st_hip) - l2*sin(q(3) +q_st_hip + q_st_knee);
    pz - l1*cos(q(3) +q_st_hip) - l2*cos(q(3) +q_st_hip + q_st_knee)
];

disp('Stance foot position:')
disp(p_st)

%% =========================
%% Swing foot position
%% =========================

p_sw = [
    px - l1*sin(q(3) +q_sw_hip) - l2*sin(q(3) +q_sw_hip + q_sw_knee);
    pz - l1*cos(q(3) +q_sw_hip) - l2*cos(q(3) +q_sw_hip + q_sw_knee)
];

disp('Swing foot position:')
disp(p_sw)

%% =========================
%% Stance Jacobian
%% =========================

J_st = simplify(jacobian(p_st,q));

disp('Stance foot Jacobian:')
disp(J_st)

%% =========================
%% Swing Jacobian
%% =========================

J_sw = simplify(jacobian(p_sw,q));

disp('Swing foot Jacobian:')
disp(J_sw)

%% =========================
%% Jdot * dq (stance)
%% =========================

Jdotdq_st = simplify(jacobian(J_st*dq,q) * dq);

disp('Jdot * dq stance:')
disp(Jdotdq_st)

%% =========================
%% Jdot * dq (swing)
%% =========================

Jdotdq_sw = simplify(jacobian(J_sw*dq,q) * dq);

disp('Jdot * dq swing:')
disp(Jdotdq_sw)

%% =========================
%% Export MATLAB functions
%% =========================

matlabFunction(J_st, ...
    'File','Dynamics/J_stance', ...
    'Vars',{q,param});

matlabFunction(Jdotdq_st, ...
    'File','Dynamics/Jdotdq_stance', ...
    'Vars',{q,dq,param});

matlabFunction(J_sw, ...
    'File','Dynamics/J_swing', ...
    'Vars',{q,param});

matlabFunction(Jdotdq_sw, ...
    'File','Dynamics/Jdotdq_swing', ...
    'Vars',{q,dq,param});

%% =========================
%% Example numerical usage
%% =========================

p = parameters();

par = packParameters(p);

q_num  = zeros(7,1);
dq_num = zeros(7,1);

J_num = J_stance(q_num,par);
Jdotdq_num = Jdotdq_stance(q_num,dq_num,par);

disp('Numerical stance Jacobian:')
disp(J_num)

disp('Numerical Jdot*dq stance:')
disp(Jdotdq_num)