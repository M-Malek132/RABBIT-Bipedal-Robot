function [value, isterminal, direction] = rabbit_impact_event(t, x)

q = x(1:7);
swing_foot = P_sw(q);
% Always return the true continuous foot height
value = swing_foot(2);

% Only stop the simulation if enough time has passed for the foot to lift off
if t < 0.001
    isterminal = 0; % Ignore zero-crossings right after step starts
else
    isterminal = 1; % Stop simulation on impact
end

% Only trigger when foot is falling
direction = -1;


end
