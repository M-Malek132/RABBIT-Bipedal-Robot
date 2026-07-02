function [stance_foot, swing_foot, hip, stance_knee, swing_knee, torso_top] = rabbit_points(q)

stance_foot = P_st(q);
swing_foot  = P_sw(q);
hip         = P_hip(q);
stance_knee = P_knee_st(q);
swing_knee  = P_knee_sw(q);
torso_top   = P_torso(q);

end
