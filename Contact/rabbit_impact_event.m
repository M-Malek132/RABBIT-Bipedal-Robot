function [value, isterminal, direction] = rabbit_impact_event(t, x, params)

q = x(1:7);
p = packParameters(params);
[~,swing_foot,~,~,~,~] = rabbit_kinematics(q, p);

% Always return the true continuous foot height
value = swing_foot(2);

% Only stop the simulation if enough time has passed for the foot to lift off
if t < 0.05
    isterminal = 0; % Ignore zero-crossings right after step starts
else
    isterminal = 1; % Stop simulation on impact
end

% Only trigger when foot is falling
direction = -1;

end
