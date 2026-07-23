function [value, isterminal, direction] = rabbit_impact_event(t, x)
%RABBIT_IMPACT_EVENT  ode45 event: swing foot strikes the ground.
%
%   [value, isterminal, direction] = rabbit_impact_event(t, x)
%
% Guard function for the hybrid single-support phase. The impact surface S is
% "swing-foot height = 0" (world frame, up-positive z), so value crossing zero
% from above marks touchdown and terminates the step.
%
%   value      = P_sw(q)_z            swing-foot height (the zero-crossing)
%   isterminal = 1 (stop) once past a short blanking window, else 0
%   direction  = -1                   only trigger while the foot is DESCENDING
%
% The t < 0.001 s blanking window suppresses a spurious immediate trigger: at
% the start of a step the swing foot is at (or numerically just below) z = 0,
% and without the guard ode45 would terminate on that initial crossing before
% the foot has lifted off. 0.001 s is long enough for lift-off but far shorter
% than a step (~1 s).

q = x(1:7);
swing_foot = P_sw(q);
% Always return the true continuous foot height (the event surface value).
value = swing_foot(2);

% Blank the first millisecond so the foot can lift off before we arm the event.
if t < 0.001
    isterminal = 0; % ignore zero-crossings right after the step starts
else
    isterminal = 1; % stop the simulation on impact
end

% Only trigger on a DOWNWARD crossing (foot descending onto the ground).
direction = -1;

end
