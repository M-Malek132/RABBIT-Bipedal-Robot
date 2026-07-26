function pts = ch3_body_points(q)
%CH3_BODY_POINTS  World-frame stick-figure points for a configuration.
%
%   pts = ch3_body_points(q)
%
% Pulls the five drawable points out of the generated homogeneous transforms.
% Column 3 of each 4x4 transform's translation is the WORLD Z (up-positive) --
% the transforms already carry the pz sign convention, so heights taken from
% here never need the -pz correction that raw generalized coordinates do.
%
% Torso length 0.75 m matches the value baked into the symbolic derivation in
% Dynamics/rabbit_energy_model_generalized_Lagrange.m. Changing it there
% without changing it here would draw a robot that is not the one simulated.
%
% Output
%   pts : struct with 2x1 fields hip, torso_top, stance_knee, swing_knee,
%         stance_foot, swing_foot   (each [x; z], world frame)
%
% See also CH3_ANIMATE, CH3_PLOT_GAIT.

TORSO_LEN = 0.75;

hip   = Tt(q) * [0;          0; 0; 1];
top   = Tt(q) * [0; -TORSO_LEN; 0; 1];
knee1 = T2(q) * [0;          0; 0; 1];
knee2 = T4(q) * [0;          0; 0; 1];

pts.hip         = [hip(1);   hip(3)];
pts.torso_top   = [top(1);   top(3)];
pts.stance_knee = [knee1(1); knee1(3)];
pts.swing_knee  = [knee2(1); knee2(3)];
pts.stance_foot = P_st(q);
pts.swing_foot  = P_sw(q);

end
