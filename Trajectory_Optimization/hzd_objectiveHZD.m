function J = hzd_objectiveHZD(z, model, opt)
%HZD_OBJECTIVEHZD  Cost: torque-squared integral / step length.

params = model.params;

[CP, q0, dq0, T] = hzd_unpackDecisionVars(z, model, opt);

simOpt.CP = CP;
simOpt.T  = T;
simOpt.Kp = opt.Kp;
simOpt.Kd = opt.Kd;

x0 = [q0; dq0];

try
    [t, x, u] = hzd_simulateOneStep(x0, model, opt, simOpt);

    if any(isnan(x(:))) || any(isnan(u(:))) || length(t) < 3
        J = 1e8; return;
    end

    nq      = model.nq;
    kin_s   = rabbit_kinematics(x(1,   1:nq)', params);
    kin_e   = rabbit_kinematics(x(end, 1:nq)', params);
    stepLen = kin_e.swingFoot(1) - kin_s.swingFoot(1);

    if stepLen <= 0.05
        J = 1e8; return;
    end

    J = trapz(t, sum(u.^2, 2)) / stepLen;

catch
    J = 1e8;
end
end
