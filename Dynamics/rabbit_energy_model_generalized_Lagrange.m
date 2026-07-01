%% +++++++++++++++++++++++++++++++++++++++++ %%
clc;clear;close all


%% +++++++++++++++++++++++++++++++++++++++++ %%


X = sym('theta',[7 1],'real');

x = X(1);   z = X(2);   qt = X(3);   q1 = X(4);  q2 = X(5);     q3 = X(6);  q4 = X(7);

% Transformation Matrixs
% forward
T_wt = [cos(qt)       -sin(qt)           0     x;
            0                   0                  1     0;
            -sin(qt)        -cos(qt)          0     -z;
            0                   0                  0     1];

T_t1 = [cos(q1)      -sin(q1)        0     0;
             sin(q1)       cos(q1)        0     0;
             0                  0                1     0;
              0                 0                0     1];

T_12 = [cos(q2)      -sin(q2)        0     0;
             sin(q2)       cos(q2)        0     0.5;
             0                  0                1     0;
            0                   0                0     1];

T_t3 = [cos(q3)      -sin(q3)        0     0;
             sin(q3)       cos(q3)        0     0;
             0                  0                1     0;
              0                 0                0     1];

T_34 = [cos(q4)      -sin(q4)        0     0;
             sin(q4)       cos(q4)        0     0.5;
             0                  0                1     0;
            0                   0                0     1];

T_tmp = [eye(3) [1 0 1.5]'; zeros(1,3) 1]* T_wt;
T_w1 = T_tmp*   T_t1;
T_w2 = T_w1*    T_12;
T_w3 = T_tmp*    T_t3;
T_w4 = T_w3*    T_34;

Pos_stance = T_w2* [0 0.5 0 1]';
Pos_swing = T_w4* [0 0.5 0 1]';

T_wt = matlabFunction(T_tmp,'File','Tt','Vars',{X});
T_w1 = matlabFunction(T_w1,'File','T1','Vars',{X});
T_w2 = matlabFunction(T_w2,'File','T2','Vars',{X});
T_w3 = matlabFunction(T_w3,'File','T3','Vars',{X});
T_w4 = matlabFunction(T_w4,'File','T4','Vars',{X});

Pos_stance = matlabFunction(Pos_stance,'File','P_st','Vars',{X});
Pos_swing = matlabFunction(Pos_swing,'File','P_sw','Vars',{X});
%% +++++++++++++++++++++++++++++++++++++++++ %%

clc;clear;close all
syms theta [7 1]

M=sym('M', [7 7]);
[~, ~, I] = Mass_Properties();

for j=1:7
    for k=1:7
         At=trace(diff(Tt(theta),theta(j))*I(:,:,1)*diff(transpose(Tt(theta)),theta(k)));
         A1=trace(diff(T1(theta),theta(j))*I(:,:,2)*diff(transpose(T1(theta)),theta(k)));
         A2=trace(diff(T2(theta),theta(j))*I(:,:,3)*diff(transpose(T2(theta)),theta(k)));
         A3=trace(diff(T3(theta),theta(j))*I(:,:,4)*diff(transpose(T3(theta)),theta(k)));
         A4=trace(diff(T4(theta),theta(j))*I(:,:,5)*diff(transpose(T4(theta)),theta(k)));
         temp = At+A1+A2+A3+A4;
         A=(temp);
         M(j,k)=A;
    end
end

M=simplify(M);

M = matlabFunction(M,'File','M','Vars',{theta});
%% +++++++++++++++++++++++++++++++++++++++++ %%

% Calculating M_dot wrt Time

clc;clear;close all
syms a(t) [7 1]
syms q [7 1]
syms Dq [7 1]

in1 = a(t);
DM=diff(M(in1),t);
DM=simplify(DM);

DM=subs(DM,[a(t) ;diff(a1, t); diff(a2, t) ;diff(a3, t); diff(a4, t); diff(a5, t); diff(a6, t); diff(a7, t)],[q; Dq] );
var = [q; Dq];
DM = matlabFunction(DM,'File','DM','var',{var});

%% +++++++++++++++++++++++++++++++++++++++++ %%

% Kinetic Energy

clc;clear
syms theta [7 1]
syms omega [7 1]

K = 1/2 * omega' * M(theta) * omega;
K = simplify(K);

% temp0  = jacobian(K,theta);
temp1 = diff(K,theta1);
temp2 = diff(K,theta2);
temp3 = diff(K,theta3);
temp4 = diff(K,theta4);
temp5 = diff(K,theta5);
temp6 = diff(K,theta6);
temp7 = diff(K,theta7);
temp  = [temp1;temp2;temp3;temp4;temp5;temp6;temp7];
in = [theta; omega];
V = (DM(in)*omega - temp);
V = simplify(V);
V = matlabFunction(V,'File','V','var',{in});

%% +++++++++++++++++++++++++++++++++++++++++ %%
% Potential Energy
clc;clear

[m, com, ~] = Mass_Properties();
syms theta [7 1]
g=[0 0 -9.8062 0]';

Pt=m(1)*g'*Tt(theta)*com(:,1);
P1=m(2)*g'*T1(theta)*com(:,2);
P2=m(3)*g'*T2(theta)*com(:,3);
P3=m(4)*g'*T3(theta)*com(:,4);
P4=m(5)*g'*T4(theta)*com(:,5);

P = (-(Pt+P1+P2+P3+P4));

G1 = diff(P,theta1);
G2 = diff(P,theta2);
G3 = diff(P,theta3);
G4 = diff(P,theta4);
G5 = diff(P,theta5);
G6 = diff(P,theta6);
G7 = diff(P,theta7);
G = [G1;G2;G3;G4;G5;G6;G7];
G = simplify(G);

G = matlabFunction(G,'File','G','var',{theta});


%% +++++++++++++++++++++++++++++++++++++++++ %%
function [m, com, I] = Mass_Properties()
% +++++++++ %% +++++++++ %% +++++++++ %

Mt = 10; %Kg

% Center of mass: ( meters )
	X = 0;
	Y = -0.75/2;
	Z = 0;
comt = [X Y Z]';

% Moments of inertia: (Kgrams *  square meters )
% Taken at the  COM
	Ixx = 0.2	;Ixy = 0	;Ixz = 0;
	Iyx = 0	;Iyy = 0.05;	Iyz = 0;
	Izx = 0	;Izy = 0	;Izz = 0.2;
    
    Ixx = Ixx + Mt*(Y^2+Z^2);
    Iyy = Iyy + Mt*(X^2+Z^2);
    Izz = Izz + Mt*(Y^2+X^2);

    I_t = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             Mt*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              Mt*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         Mt*Z; ...
            Mt*X                Mt*Y                Mt*Z             Mt];

    % +++++++++ %% +++++++++ %% +++++++++ %


M1 = 5; %Kg

% Center of mass: ( meters )
	X = 0;
	Y = 0.5/2;
	Z = 0;
com1 = [X Y Z]';

% Moments of inertia: (Kgrams *  square meters )
% Taken at the COM
	Ixx = 0.1	;Ixy = 0	;Ixz = 0;
	Iyx = 0	;Iyy = 0.02;	Iyz = 0;
	Izx = 0	;Izy = 0	;Izz = 0.1;
    
    Ixx = Ixx + M1*(Y^2+Z^2);
    Iyy = Iyy + M1*(X^2+Z^2);
    Izz = Izz + M1*(Y^2+X^2);

    I_1 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M1*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M1*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M1*Z; ...
            M1*X                M1*Y                M1*Z             M1];

% +++++++++ %% +++++++++ %% +++++++++ %

M2 = 5; %Kg

% Center of mass: (meters )
	X = 0;
	Y = 0.5/2;
	Z = 0;
com2 = [X Y Z]';

% Moments of inertia: (Kgrams *  square meters )
% Taken at the COM
	Ixx = 0.1	;Ixy = 0	;Ixz = 0;
	Iyx = 0	;Iyy = 0.02;	Iyz = 0;
	Izx = 0	;Izy = 0	;Izz = 0.1;
    
    Ixx = Ixx + M2*(Y^2+Z^2);
    Iyy = Iyy + M2*(X^2+Z^2);
    Izz = Izz + M2*(Y^2+X^2);

    I_2 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M2*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M2*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M2*Z; ...
            M2*X                M2*Y                M2*Z             M2];

    % +++++++++ %% +++++++++ %% +++++++++ %


M3 = 5; %Kg

% Center of mass: ( meters )
	X = 0;
	Y = 0.5/2;
	Z = 0;
com3 = [X Y Z]';

% Moments of inertia: (Kgrams *  square meters )
% Taken at the COM
	Ixx = 0.1	;Ixy = 0	;Ixz = 0;
	Iyx = 0	;Iyy = 0.02;	Iyz = 0;
	Izx = 0	;Izy = 0	;Izz = 0.1;
    
    Ixx = Ixx + M3*(Y^2+Z^2);
    Iyy = Iyy + M3*(X^2+Z^2);
    Izz = Izz + M3*(Y^2+X^2);

    I_3 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M3*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M3*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M3*Z; ...
            M3*X                M3*Y                M3*Z             M3];

% +++++++++ %% +++++++++ %% +++++++++ %

M4 = 5; %Kg

% Center of mass: (meters )
	X = 0;
	Y = 0.5/2;
	Z = 0;
com4 = [X Y Z]';

% Moments of inertia: (Kgrams *  square meters )
% Taken at the COM
	Ixx = 0.1	;Ixy = 0	;Ixz = 0;
	Iyx = 0	;Iyy = 0.02;	Iyz = 0;
	Izx = 0	;Izy = 0	;Izz = 0.1;
    
    Ixx = Ixx + M4*(Y^2+Z^2);
    Iyy = Iyy + M4*(X^2+Z^2);
    Izz = Izz + M4*(Y^2+X^2);

    I_4 = [	(-Ixx+Iyy+Izz)/2    Ixy                 Ixz             M4*X; ...
	        Iyx             (Ixx-Iyy+Izz)/2        Iyz              M4*Y; ...
	        Izx                 Izy         (Ixx+Iyy-Izz)/2         M4*Z; ...
            M4*X                M4*Y                M4*Z             M4];

m = [Mt M1 M2 M3 M4]';
com = [comt com1 com2 com3 com4];
com = [com; ones(1,5)];
I(:,:,1) = I_t;     I(:,:,2) = I_1;     I(:,:,3) = I_2;     I(:,:,4) = I_3;     I(:,:,5) = I_4;

end
