function [cineq, ceq] = hzd_constraintsHZD(z, model, opt)
%HZD_CONSTRAINTSHZD  Nonlinear constraints for B-spline HZD optimisation.
%
%  Equality (ceq = 0):
%    1. Periodicity after impact + relabelling
%    2. Desired average walking speed
%    3. Swing foot height = 0 at impact
%    4. HZD invariance: y=0, dy=0 at start; y=0, dy=0 after impact
%
%  Inequality (cineq <= 0):
%    Swing-foot clearance, hip height, torque limits, joint limits, knees

params = model.params;
nq     = model.nq;

[CP, q0, dq0, T] = hzd_unpackDecisionVars(z, model, opt);

simOpt.CP = CP;
simOpt.T  = T;
simOpt.Kp = opt.Kp;
simOpt.Kd = opt.Kd;

cineq = [];
ceq   = [];

try
    [t, x, u] = hzd_simulateOneStep([q0; dq0], model, opt, simOpt);

    if any(isnan(x(:))) || any(isnan(u(:))) || length(t) < 3
        cineq = 1e6*ones(200,1); ceq = 1e6*ones(30,1); return;
    end

    qStart  = x(1,   1:nq)';   dqStart = x(1,   nq+1:end)';
    qEnd    = x(end, 1:nq)';   dqEnd   = x(end, nq+1:end)';

    %% 1. Impact + relabelling
    xPlus      = rabbit_impact_map([qEnd; dqEnd], params);
    xRelabeled = rabbit_reset_map(xPlus, params);
    qRel  = xRelabeled(1:nq);
    dqRel = xRelabeled(nq+1:end);

    %% 2. Periodicity
    ceq = [ceq; qRel - qStart; dqRel - dqStart];

    %% 3. Speed
    kin_s   = rabbit_kinematics(qStart, params);
    kin_e   = rabbit_kinematics(qEnd,   params);
    stepLen = kin_e.swingFoot(1) - kin_s.swingFoot(1);
    ceq = [ceq; stepLen/T - opt.v_des];

    %% 4. Foot height at impact
    ceq = [ceq; kin_e.swingFoot(2)];

    %% 5. HZD invariance
    [y0,    dy0   ] = hzd_virtualConstraints(qStart, dqStart, CP, model, opt);
    [yPlus, dyPlus] = hzd_virtualConstraints(qRel,   dqRel,  CP, model, opt);
    ceq = [ceq; y0; dy0; yPlus; dyPlus];

    %% 6. Path inequality constraints
    N = length(t);
    footClr   = zeros(N,1);
    hipViol   = zeros(N,1);
    torqueHi  = zeros(N, model.nu);
    torqueLo  = zeros(N, model.nu);
    jointHi   = zeros(N, nq);
    jointLo   = zeros(N, nq);
    kneeViol  = zeros(N,2);

    for k = 1:N
        qk  = x(k, 1:nq)';
        uk  = u(k,:)';
        kin = rabbit_kinematics(qk, params);

        footClr(k)     = -kin.swingFoot(2);
        hipViol(k)     = opt.hipHeightMin - kin.hip(2);
        jointHi(k,:)   = qk' - opt.qMax';
        jointLo(k,:)   = opt.qMin' - qk';
        kneeViol(k,1)  = qk(5);   % stance knee <= 0 when bent
        kneeViol(k,2)  = qk(7);   % swing  knee <= 0 when bent
        torqueHi(k,:)  = uk' - opt.uMax;
        torqueLo(k,:)  = opt.uMin - uk';
    end

    cineq = [footClr; hipViol; jointHi(:); jointLo(:); ...
             kneeViol(:); torqueHi(:); torqueLo(:)];

catch ME
    warning('hzd_constraintsHZD: %s', ME.message);
    cineq = 1e6*ones(200,1);
    ceq   = 1e6*ones(30,1);
end
end
