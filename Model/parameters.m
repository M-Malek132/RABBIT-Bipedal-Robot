function p = parameters()
% PARAMETERS  Return structure of physical parameters for RABBIT robot.

% Masses [kg]
p.mT = 10;     % torso
p.m1 = 5;      % thigh
p.m2 = 5;      % shank

% Link lengths [m]
p.l1 = 0.5;
p.l2 = 0.5;
p.lt = 0.75;

% Inertias [kg·m²]
p.I1 = 0.1;
p.I2 = 0.1;
p.IT = 0.2;

% Gravity [m/s²]
p.g  = 9.81;
end
