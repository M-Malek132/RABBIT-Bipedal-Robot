function [t, x, uHist] = hzd_simulateOneStep(x0, model, opt, simOpt)
%HZD_SIMULATEONESTEP  Simulate one swing phase using HZD controller.
%
%  Uses hzd_closedLoopDynamics and hzd_footStrikeEvent.
%  Does NOT call your repo's simulate_one_step.m directly because
%  that function may use a different controller interface.
%  If you prefer to reuse simulate_one_step.m, wrap it here instead.

T = simOpt.T;

odeOpts = odeset( ...
    'RelTol',    1e-6, ...
    'AbsTol',    1e-8, ...
    'Events',    @(t,x) hzd_footStrikeEvent(t, x, model, opt), ...
    'MaxStep',   0.01);

[t, x] = ode45( ...
    @(t,x) hzd_closedLoopDynamics(t, x, model, opt, simOpt), ...
    [0, T], x0, odeOpts);

% Reconstruct torque history
n     = length(t);
uHist = zeros(n, model.nu);
for k = 1:n
    [~, uk] = hzd_closedLoopDynamics(t(k), x(k,:)', model, opt, simOpt);
    uHist(k,:) = uk';
end
end
