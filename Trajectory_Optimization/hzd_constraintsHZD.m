function [cineq, ceq] = hzd_constraintsHZD(z, model, opt)
%HZD_CONSTRAINTSHZD Nonlinear constraints for B-spline HZD optimisation.
%
% Equality constraints ceq = 0:
%   1. Phase starts and ends at desired values
%   2. Periodicity after impact + relabelling
%   3. Desired average walking speed
%   4. Swing foot height = 0 at impact
%   5. HZD invariance: y = 0, dy = 0 at start and after impact
%
% Inequality constraints cineq <= 0:
%   Swing-foot clearance, hip height, torque limits, joint limits, knees

    params = model.params;
    p      = packParameters(params);

    nq = model.nq;
    nu = model.nu;

    penalty = 1e6;

    cineq = [];
    ceq   = [];

    try
        %% Unpack decision variables
        [CP, q0, dq0, T] = hzd_unpackDecisionVars(z, model, opt);

        q0  = q0(:);
        dq0 = dq0(:);

        if numel(q0) ~= nq
            error('q0 has wrong number of elements. Expected %d, got %d.', nq, numel(q0));
        end

        if numel(dq0) ~= nq
            error('dq0 has wrong number of elements. Expected %d, got %d.', nq, numel(dq0));
        end

        x0 = [q0; dq0];

        %% Phase start equality
        [theta0, ~] = hzd_phaseVariable(q0, model);
        ceq = [ceq; theta0 - opt.thetaStart];

        %% Simulation options
        simOpt.CP = CP;
        simOpt.T  = T;
        simOpt.Kp = opt.Kp;
        simOpt.Kd = opt.Kd;

        %% Simulate one step
        [t, x, u] = hzd_simulateOneStep(x0, model, opt, simOpt);
        animate_rabbit_stepping_stones(x,params);
        %% Validate and normalize simulation output
        [t, x, u, isValid] = normalizeSimulationOutput(t, x, u, nq, nu);

        if ~isValid
            [cineq, ceq] = failedConstraints(penalty);
            return;
        end

        if any(isnan(x(:))) || any(isnan(u(:))) || numel(t) < 3
            [cineq, ceq] = failedConstraints(penalty);
            return;
        end

        %% Extract start/end states
        qStart  = x(1,   1:nq).';
        dqStart = x(1,   nq+1:2*nq).';

        qEnd    = x(end, 1:nq).';
        dqEnd   = x(end, nq+1:2*nq).';

        %% Phase end equality
        [thetaEnd, ~] = hzd_phaseVariable(qEnd, model);
        ceq = [ceq; thetaEnd - opt.thetaEnd];

        %% Impact + relabelling
        xRelabeled = rabbit_reset_map([qEnd; dqEnd], params);

        qRel  = xRelabeled(1:nq);
        dqRel = xRelabeled(nq+1:2*nq);

        %% Periodicity
        ceq = [ceq;
               qRel  - qStart;
               dqRel - dqStart];

        %% Desired speed
        [~, swingFootStart, ~, ~, ~, ~] = rabbit_kinematics(qStart, p);
        [~, swingFootEnd,   ~, ~, ~, ~] = rabbit_kinematics(qEnd,   p);

        stepLen = swingFootEnd(1) - swingFootStart(1);

        ceq = [ceq;
               stepLen / T - opt.v_des];

        %% Swing foot height at impact
        ceq = [ceq;
               swingFootEnd(2)];

        %% HZD invariance
        [y0,    dy0]    = hzd_virtualConstraints(qStart, dqStart, CP, model, opt);
        [yPlus, dyPlus] = hzd_virtualConstraints(qRel,   dqRel,   CP, model, opt);

        ceq = [ceq;
               y0;
               dy0;
               yPlus;
               dyPlus];

        %% Path inequality constraints
        N = size(x, 1);

        footClr  = zeros(N, 1);
        hipViol  = zeros(N, 1);
        jointHi  = zeros(N, nq);
        jointLo  = zeros(N, nq);
        kneeViol = zeros(N, 2);
        torqueHi = zeros(N, nu);
        torqueLo = zeros(N, nu);

        qMax = opt.qMax(:).';
        qMin = opt.qMin(:).';
        uMax = opt.uMax(:).';
        uMin = opt.uMin(:).';

        for k = 1:N
            qk = x(k, 1:nq).';
            uk = u(k, :).';
            
            [~, swingFoot, hip, ~, ~, ~] = rabbit_kinematics(qk, p);

            footClr(k)    = -swingFoot(2);
            hipViol(k)    = opt.hipHeightMin - hip(2);

            jointHi(k, :) = qk.' - qMax;
            jointLo(k, :) = qMin - qk.';

            kneeViol(k, 1) = qk(5);
            kneeViol(k, 2) = qk(7);

            torqueHi(k, :) = uk.' - uMax;
            torqueLo(k, :) = uMin - uk.';
        end

        cineq = [footClr;
                 hipViol;
                 jointHi(:);
                 jointLo(:);
                 kneeViol(:);
                 torqueHi(:);
                 torqueLo(:)];

    catch ME
        fprintf(2, '\n[hzd_constraintsHZD ERROR]\n');
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        
        [cineq, ceq] = failedConstraints(penalty);
    end
end


function [t, x, u, isValid] = normalizeSimulationOutput(t, x, u, nq, nu)
%NORMALIZESIMULATIONOUTPUT Make simulation outputs consistent.
%
% Expected:
%   t : N x 1
%   x : N x 2*nq
%   u : N x nu

    isValid = false;

    if isempty(t) || isempty(x) || isempty(u)
        return;
    end

    t = t(:);

    % Normalize x to N x 2*nq
    if size(x, 2) ~= 2*nq && size(x, 1) == 2*nq
        x = x.';
    end

    % Normalize u to N x nu
    if size(u, 2) ~= nu && size(u, 1) == nu
        u = u.';
    end

    if size(x, 2) ~= 2*nq
        return;
    end

    if size(u, 2) ~= nu
        return;
    end

    if size(x, 1) ~= size(u, 1)
        return;
    end

    if numel(t) ~= size(x, 1)
        return;
    end

    isValid = true;
end


function [cineq, ceq] = failedConstraints(penalty)
%FAILEDCONSTRAINTS Return large constraint violations after simulation failure.

    cineq = penalty * ones(200, 1);
    ceq   = penalty * ones(30, 1);
end
