function [value,isterminal,direction] = impact_event_wrapper(t,xi)
    nq = 7;
    [value,isterminal,direction] = rabbit_impact_event(t, xi(1:2*nq));
end
