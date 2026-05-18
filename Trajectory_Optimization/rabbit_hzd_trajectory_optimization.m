%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RABBIT HZD TRAJECTORY OPTIMIZATION
%
% This example demonstrates:
%
% 1. Direct collocation formulation
% 2. Bezier trajectory parameterization
% 3. Cost function definition
% 4. Gait periodicity constraints
% 5. Nonlinear optimization using fmincon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%% =========================================================
% PARAMETERS
%% =========================================================

N  = 30;          % collocation nodes
nx = 10;          % state dimension [q(5); dq(5)]
nu = 4;           % actuators
dt = 0.03;        % timestep

p.N  = N;
p.nx = nx;
p.nu = nu;
p.dt = dt;

%% =========================================================
% INITIAL GUESS
%% =========================================================

X0 = zeros(nx,N);      % states
U0 = zeros(nu,N-1);    % controls

z0 = [X0(:); U0(:)];

%% =========================================================
% OPTIMIZATION
%% =========================================================

options = optimoptions('fmincon',...
    'Display','iter',...
    'MaxFunctionEvaluations',1e6,...
    'MaxIterations',500,...
    'Algorithm','sqp');

[z_opt,fval] = fmincon(...
    @(z) cost_function(z,p),...
    z0,...
    [],[],[],[],[],[],...
    @(z) gait_constraints(z,p),...
    options);

%% =========================================================
% EXTRACT SOLUTION
%% =========================================================

[X,U] = unpack_decision_variables(z_opt,p);

disp('Optimization complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COST FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = cost_function(z,p)

[X,U] = unpack_decision_variables(z,p);

J = 0;

for k = 1:p.N-1

    u = U(:,k);

    % torque squared cost
    J = J + u'*u;

end

J = J * p.dt;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAIT CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,ceq] = gait_constraints(z,p)

[X,U] = unpack_decision_variables(z,p);

ceq = [];
c   = [];

%% =========================================================
% DIRECT COLLOCATION CONSTRAINTS
%% =========================================================

for k = 1:p.N-1

    xk   = X(:,k);
    xkp1 = X(:,k+1);

    uk = U(:,k);

    % dynamics at node k
    fk = rabbit_dynamics(xk,uk);

    % Euler collocation
    x_next = xk + p.dt * fk;

    ceq = [ceq;
           xkp1 - x_next];

end

%% =========================================================
% PERIODICITY CONSTRAINT
%% =========================================================

x0 = X(:,1);
xF = X(:,end);

% impact map
x_plus = rabbit_impact_map(xF);

% periodic gait
ceq = [ceq;
       x0 - x_plus];

%% =========================================================
% FOOT CONTACT CONSTRAINTS
%% =========================================================

for k = 1:p.N

    q = X(1:5,k);

    [~,swingFoot] = foot_positions(q);

    % swing foot height >= 0
    c = [c;
         -swingFoot(2)];

end

%% =========================================================
% IMPACT CONSTRAINT
%% =========================================================

qF = X(1:5,end);

[~,swingFoot] = foot_positions(qF);

% swing foot touches ground
ceq = [ceq;
       swingFoot(2)];

%% =========================================================
% TORQUE LIMITS
%% =========================================================

umax = 50;

for k = 1:p.N-1

    u = U(:,k);

    c = [c;
         u - umax;
        -u - umax];

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = rabbit_dynamics(x,u)

% states
q  = x(1:5);
dq = x(6:10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% robot dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = rabbit_D(q);
C = rabbit_C(q,dq);
G = rabbit_G(q);
B = rabbit_B();

ddq = D \ (B*u - C*dq - G);

dx = [dq;
      ddq];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNPACK VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,U] = unpack_decision_variables(z,p)

nx = p.nx;
nu = p.nu;
N  = p.N;

Xsize = nx*N;

X = reshape(z(1:Xsize),nx,N);

U = reshape(z(Xsize+1:end),nu,N-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPACT MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_plus = rabbit_impact_map(x_minus)

q  = x_minus(1:5);
dq = x_minus(6:10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simplified impact model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dq_plus = 0.8 * dq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% leg relabeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_plus = relabel(q);

x_plus = [q_plus;
          dq_plus];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RELABELING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q_new = relabel(q)

q_new = q;

% swap legs example
tmp      = q_new(2);
q_new(2) = q_new(4);
q_new(4) = tmp;

end
