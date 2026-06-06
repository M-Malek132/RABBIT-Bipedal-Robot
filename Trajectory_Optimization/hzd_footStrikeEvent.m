function [value, isterminal, direction] = hzd_footStrikeEvent(t, x, model, opt)
%HZD_FOOTSTRIKEEVENT  ODE event: swing foot touches ground.
%
%  Calls rabbit_kinematics(q, params) from your Model/ folder.
%  Expected output: kin.swingFoot = [x; z]  (position of swing foot tip)
%
%  If your function uses different field names, change kin.swingFoot below.

nq     = model.nq;
params = model.params;

q = x(1:nq);

kin = rabbit_kinematics(q, params);

% Swing foot vertical height — should cross zero at touchdown
value      = kin.swingFoot(2);

isterminal = 1;     % stop integration
direction  = -1;    % detect only downward crossing

% Suppress false detection at t = 0
if t < 0.04
    value = 1;
end
end
