% ============================================================
% File:
% dynamics/input_matrix.m
% ============================================================

function B = input_matrix()

% ------------------------------------------------------------
% Generalized coordinates:
%
% q = [x z qt q1 q2 q3 q4]'
%
% x   : torso horizontal position
% z   : torso vertical position
% qt  : torso pitch angle
% q1  : stance hip
% q2  : stance knee
% q3  : swing hip
% q4  : swing knee
%
% Total DOF = 7
%
% ------------------------------------------------------------
% RABBIT has 4 actuators:
%
% u = [u1 u2 u3 u4]'
%
% Actuated joints:
%   q1 q2 q3 q4
%
% Unactuated coordinates:
%   x z qt
%
% ------------------------------------------------------------
% Dynamics:
%
% D(q)qdd + C(q,dq)dq + G(q) = B*u
%
% ------------------------------------------------------------

B = [ 0  0  0  0 ;
      0  0  0  0 ;
      0  0  0  0 ;
      1  0  0  0 ;
      0  1  0  0 ;
      0  0  1  0 ;
      0  0  0  1 ];

end
