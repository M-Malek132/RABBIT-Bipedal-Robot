%% +++++++++++++++++++++++++++++++++++++++++ %%
clc;clear;close all


%% +++++++++++++++++++++++++++++++++++++++++ %%


X = sym('theta',[3 1],'real');

x = X(1);   z = X(2);   qt = X(3);  

% Transformation Matrixs
% forward
T_wt = [cos(qt)       -sin(qt)           0     x;
            0                   0                  1     0;
            -sin(qt)        -cos(qt)          0     -z;
            0                   0                  0     1];


T_tmp = [eye(3) [1 0 1.5]'; zeros(1,3) 1]* T_wt;


T_wt = matlabFunction(T_tmp,'File','Tt','Vars',{X});


%% +++++++++++++++++++++++++++++++++++++++++ %%

clc;clear;close all
syms theta [3 1]

M=sym('M', [3 3]);
[~, ~, I] = Mass_Properties();

for j=1:3
    for k=1:3
         At=trace(diff(Tt(theta),theta(j))*I(:,:,1)*diff(transpose(Tt(theta)),theta(k)));
         temp = At;
         A=(temp);
         M(j,k)=A;
    end
end

M=simplify(M);

M = matlabFunction(M,'File','M','Vars',{theta});

%% +++++++++++++++++++++++++++++++++++++++++ %%

% Calculating M_dot wrt Time

clc;clear;close all
syms a(t) [3 1]
syms q [3 1]
syms Dq [3 1]

in1 = a(t);
DM=diff(M(in1),t);
DM=simplify(DM);

DM=subs(DM,[a(t) ;diff(a1, t); diff(a2, t) ;diff(a3, t)],[q; Dq] );
var = [q; Dq];
DM = matlabFunction(DM,'File','DM','var',{var});

%% +++++++++++++++++++++++++++++++++++++++++ %%

% Kinetic Energy

clc;clear
syms theta [3 1]
syms omega [3 1]

K = 1/2 * omega' * M(theta) * omega;
K = simplify(K);

% temp0  = jacobian(K,theta);
temp1 = diff(K,theta1);
temp2 = diff(K,theta2);
temp3 = diff(K,theta3);
temp  = [temp1;temp2;temp3];
in = [theta; omega];
V = (DM(in)*omega - temp);
V = simplify(V);
V = matlabFunction(V,'File','V','var',{in});

%% +++++++++++++++++++++++++++++++++++++++++ %%
% Potential Energy
clc;clear

[m, com, ~] = Mass_Properties();
syms theta [3 1]
g=[0 0 -9.8062 0]';

Pt=m(1)*g'*Tt(theta)*com(:,1);

P = (-(Pt));

G1 = diff(P,theta1);
G2 = diff(P,theta2);
G3 = diff(P,theta3);
G = [G1;G2;G3];
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


m = [Mt ]';
com = [comt ];
com = [com; ones(1,1)];
I(:,:,1) = I_t;   

end
