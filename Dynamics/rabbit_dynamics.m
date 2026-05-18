function ddq = rabbit_dynamics(q,dq,u,param)
%RABBIT_DYNAMICS
% Computes unconstrained forward dynamics of the RABBIT robot
%
% Equation of motion:
%
%   D(q)ddq + C(q,dq) + G(q) = B*u
%
% Inputs:
%   q     : generalized coordinates      [7x1]
%   dq    : generalized velocities       [7x1]
%   u     : actuator torques             [4x1]
%   param : robot parameter structure/vector
%
% Output:
%   ddq   : generalized accelerations    [7x1]

%% Dynamics matrices

D = D_matrix(q,param);

C = C_vector(q,dq,param);

G = G_vector(q,param);

B = input_matrix();

%% Forward dynamics

ddq = D \ (B*u - C - G);

end
