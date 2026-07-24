%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUND VALIDITY UTILITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pts = get_body_points(q)
    pos_hip   = Tt(q) * [0 0 0 1]';
    P_torso   = Tt(q) * [0 -0.75 0 1]';
    P_knee_st = T2(q) * [0 0 0 1]';
    P_knee_sw = T4(q) * [0 0 0 1]';

    pts.hip         = [pos_hip(1,1);   pos_hip(3,1)];
    pts.torso_top   = [P_torso(1,1);   P_torso(3,1)];
    pts.stance_knee = [P_knee_st(1,1); P_knee_st(3,1)];
    pts.swing_knee  = [P_knee_sw(1,1); P_knee_sw(3,1)];
    pts.stance_foot = P_st(q);
    pts.swing_foot  = P_sw(q);
end