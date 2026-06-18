% 
%     Mass properties of L1
%      Configuration: Default
%      Coordinate system: L1FCS

M1 = 16691.00*10^-3 %grams

% Center of mass: ( millimeters )
	X = -60.19*10^-3
	Y = -53.63*10^-3
	Z = -20.07*10^-3
com1 = [X Y Z]';

% Moments of inertia: ( grams *  square millimeters )
% Taken at the output coordinate system.
	Ixx = 399228807.86*10^-9	;Ixy = 54160086.30*10^-9	;Ixz = 22051154.53*10^-9;
	Iyx = 54160086.30*10^-9	;Iyy = 419257715.82*10^-9;	Iyz = 23674125.68*10^-9;
	Izx = 22051154.53*10^-9	;Izy = 23674125.68*10^-9	;Izz = 160243767.89*10^-9;
    
    I_1 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M1*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M1*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M1*Z; ...
            M1*X                M1*Y                M1*Z             M1]
    

% Mass properties of L2
%      Configuration: Default
%      Coordinate system: L2FCS

M2 = 4.50500000 %kilograms

% Center of Mass (user-overridden): ( meters )
	X = -0.00308490
	Y = -0.04997900
	Z = -0.11318900
com2 = [X Y Z]';

% Moments of inertia: ( kilograms * square meters )
% Taken at the output coordinate system.
	Ixx = 0.09328452	;Ixy =0.00069118	;Ixz = 0.00157560;
	Iyx = 0.00069118	;Iyy = 0.07427853	;Iyz = 0.03008645;
	Izx = 0.00157560	;Izy = 0.03008645	;Izz = 0.02789887;
    
    I_2 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M2*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M2*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M2*Z; ...
            M2*X                M2*Y                M2*Z             M2]
    

% Mass properties of L3
%      Configuration: Default
%      Coordinate system: L3FCS

M3 = 1.72510000 %kilograms

% Center of Mass (user-overridden): ( meters )
	X = -0.00016232
	Y = 0.02590400
	Z = -0.12555000
com3 = [X Y Z]';
% Moments of inertia: ( kilograms * square meters )
% Taken at the output coordinate system.
	Ixx = 0.04216261	;Ixy = -0.00000726	;Ixz =  0.00003516;
	Iyx = -0.00000726	;Iyy = 0.04126381	;Iyz = -0.00610850;
	Izx = 0.00003516	;Izy = -0.00610850	;Izz = 0.00322705;
    
    I_3 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M3*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M3*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M3*Z; ...
            M3*X                M3*Y                M3*Z             M3]
    

% Mass properties of L4
%      Configuration: Default
%      Coordinate system: L4FCS

M4 = 14332.00*10^-3 %grams

% Center of mass: ( millimeters )
	X = -305.44*10^-3
	Y = -87.04*10^-3
	Z = -0.20*10^-3
com4 = [X Y Z]';
% Moments of inertia: ( grams *  square millimeters )
% Taken at the output coordinate system.
	Ixx = 167942021.86*10^-9	;Ixy = 452656368.65*10^-9	;Ixz = 193465.05*10^-9;
	Iyx = 452656368.65*10^-9	;Iyy = 1824591497.45*10^-9	;Iyz = -60634.71*10^-9;
	Izx = 193465.05*10^-9	;Izy = -60634.71*10^-9	;Izz = 1977494768.42*10^-9;
    
    I_4 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M4*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M4*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M4*Z; ...
            M4*X                M4*Y                M4*Z             M4]
    

m = [M1 M2 M3 M4]';
com = [com1 com2 com3 com4];
com = [com; ones(1,4)];
I(:,:,1) = I_1;
I(:,:,2) = I_2;
I(:,:,3) = I_3;
I(:,:,4) = I_4;

clc;clear;

X = sym('theta',[4 1],'real')
x = X/2;

s1 = eye(4);
i1 = [0,-1,0,0;1,0,0,0;0,0,0,1;0,0,-1,0];
j1 = [0,0,-1,0;0,0,0,-1;1,0,0,0;0,1,0,0];
k1 = [0,0,0,-1;0,0,1,0;0,-1,0,0;1,0,0,0];

s2 = eye(4);
i2 = [0,-1,0,0;1,0,0,0;0,0,0,-1;0,0,1,0];
j2 = [0,0,-1,0;0,0,0,1;1,0,0,0;0,-1,0,0];
k2 = [0,0,0,-1;0,0,-1,0;0,1,0,0;1,0,0,0];
% Links lenghts
d0 = [0 0 80]*10^-3;
d1 = [81 50 236]*10^-3;
d2 = [0 141.5 140]*10^-3;
d3 = [0 -41.25 259]*10^-3;
d4 = [532.5 9.80 0]*10^-3;

% Frames in quaternion forms
% world to base
q0a = cos(0)*s1+sin(0)*k1;
q0b = cos(0)*s2+sin(0)*k2;
R0 = (q0b'*q0a)';
q0 = [cos(0),sin(0)*[0 0 1]]'
p0 = R0*[0,d0]'
% base to link 1
q1a = cos(x(1))*s1+sin(x(1))*k1;
q1b = cos(x(1))*s2+sin(x(1))*k2;
R1 = simplify(q1b'*q1a)';
q1 = [cos(x(1)),sin(x(1))*[0  0 1]]'
p1 = R1*[0,d1]'
% link 1 to link 2
q2a = cos(x(2))*s1+sin(x(2))*j1;
q2b = cos(x(2))*s2+sin(x(2))*j2;
R2 = simplify(q2b'*q2a)';
q2 = [cos(x(3)),sin(x(3))*[0 1 0]]'
p2 = R2*[0,d2]'
% link 2 to link 3
q3a = cos(x(3))*s1+sin(x(3))*k1;
q3b = cos(x(3))*s2+sin(x(3))*k2;
R3 = simplify(q3b'*q3a)';
q3 = [cos(x(4)),sin(x(4))*[0 0 1]]'
p3 = R3*[0,d3]'
% link 3 to link4 (End Effector)
q4a = cos(x(4))*s1-sin(x(4))*j1;
q4b = cos(x(4))*s2-sin(x(4))*j2;
R4 = simplify(q4b'*q4a)';
q4 = [cos(x(4)),sin(x(4))*[0 -1 0]]'
p4 = R4*[0,d4]'

% Transformation Matrixs
% forward
T_w0 = [R0,p0;0,0,0,0,1]
T_01 = [R1,p1;0,0,0,0,1]
T_12 = [R2,p2;0,0,0,0,1]
T_23 = [R3,p3;0,0,0,0,1]
T_34 = [R4,p4;0,0,0,0,1]

R_q = simplify(q0a*q1a*q2a*q3a*q4);

T_w1 = T_w0*T_01;
T_w2 = T_w1*T_12;
T_w3 = T_w2*T_23;
T_w4 = T_w3*T_34;

T_w1 = matlabFunction(T_w1(2:end,2:end),'File','T1','Vars',{X});
T_w2 = matlabFunction(T_w2(2:end,2:end),'File','T2','Vars',{X});
T_w3 = matlabFunction(T_w3(2:end,2:end),'File','T3','Vars',{X});
T_w4 = matlabFunction(T_w4(2:end,2:end),'File','T4','Vars',{X});

clc;clear;close all
syms theta [4 1]

M=sym('M', [4 4]);
run('massProperties2.mlx');

for j=1:4
    for k=1:4
         A1=trace(diff(T1(theta),theta(j))*I(:,:,1)*diff(transpose(T1(theta)),theta(k)));
         A2=trace(diff(T2(theta),theta(j))*I(:,:,2)*diff(transpose(T2(theta)),theta(k)));
         A3=trace(diff(T3(theta),theta(j))*I(:,:,3)*diff(transpose(T3(theta)),theta(k)));
         A4=trace(diff(T4(theta),theta(j))*I(:,:,4)*diff(transpose(T4(theta)),theta(k)));
         temp=A1+A2+A3+A4;
         A=(temp);
         M(j,k)=A;
    end
end

% M=simplify(M);

M = matlabFunction(M,'File','M','Vars',{theta})


% Calculating M_dot wrt Time

clc;clear;close all
syms a(t) [4 1]
syms q [4 1]
syms Dq [4 1]
syms m [4 1]
in1 = [a(t); m]
DM=diff(M(in1),t);
DM=(DM);

DM=subs(DM,[a(t) ;diff(a1, t); diff(a2, t) ;diff(a3, t); diff(a4, t)],[q; Dq] );
var = [q; Dq]
DM = matlabFunction(DM,'File','DM','var',{var})

% Kinetic Energy

clc;clear
syms theta [4 1]
syms omega [4 1]

K = 1/2 * omega' * M(theta) * omega;
% K = simplify(K)

% temp0  = jacobian(K,theta);
temp1 = diff(K,theta1);
temp2 = diff(K,theta2);
temp3 = diff(K,theta3);
temp4 = diff(K,theta4);
temp  = [temp1;temp2;temp3;temp4];
in = [theta; omega]
V = (DM(in)*omega - temp);
V = matlabFunction(V,'File','V','var',{in});

% Potential Energy
clc;clear

% run("massProperties2.mlx")
% syms theta [4 1]
g=[0 0 -9.8062 0]';

P1=m(1)*g'*T1(theta)*com(:,1);
P2=m(2)*g'*T2(theta)*com(:,2);
P3=m(3)*g'*T3(theta)*com(:,3);
P4=m(4)*g'*T4(theta)*com(:,4);

P = (-(P1+P2+P3+P4));
G1 = diff(P,theta1);
G2 = diff(P,theta2);
G3 = diff(P,theta3);
G4 = diff(P,theta4);
G = [G1;G2;G3;G4];

G = matlabFunction(G,'File','G','var',{theta})


    

