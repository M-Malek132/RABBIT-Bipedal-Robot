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

fprintf('\n==============================\n');
fprintf('ENTERING hzd_constraintsHZD\n');
fprintf('==============================\n');

fprintf('size(z) = [%d %d]\n', size(z,1), size(z,2));
fprintf('model.nq = %d\n', model.nq);
fprintf('model.nu = %d\n', model.nu);

params = model.params;
p = packParameters(params);
nq     = model.nq;
nu     = model.nu;

[CP, q0, dq0, T] = hzd_unpackDecisionVars(z, model, opt);

fprintf('\n--- After unpacking z ---\n');
fprintf('size(CP)  = [%d %d]\n', size(CP,1), size(CP,2));
fprintf('size(q0)  = [%d %d]\n', size(q0,1), size(q0,2));
fprintf('size(dq0) = [%d %d]\n', size(dq0,1), size(dq0,2));
fprintf('size(T)   = [%d %d]\n', size(T,1), size(T,2));
fprintf('T = %.8f\n', T);

if numel(q0) ~= model.nq
    error('q0 has wrong number of elements. Expected %d, got %d.', model.nq, numel(q0));
end

if numel(dq0) ~= model.nq
    error('dq0 has wrong number of elements. Expected %d, got %d.', model.nq, numel(dq0));
end

q0 = q0(:);
dq0 = dq0(:);
x0 = [q0; dq0];

fprintf('size(x0) = [%d %d]\n', size(x0,1), size(x0,2));

simOpt.CP = CP;
simOpt.T  = T;
simOpt.Kp = opt.Kp;
simOpt.Kd = opt.Kd;

cineq = [];
ceq   = [];

try
    fprintf('\n--- Calling hzd_simulateOneStep ---\n');

    [t, x, u] = hzd_simulateOneStep(x0, model, opt, simOpt);

    fprintf('\n--- Raw simulation outputs ---\n');
    fprintf('size(t) = [%d %d]\n', size(t,1), size(t,2));
    fprintf('size(x) = [%d %d]\n', size(x,1), size(x,2));
    fprintf('size(u) = [%d %d]\n', size(u,1), size(u,2));

    if isempty(t)
        error('Simulator returned empty t.');
    end

    if isempty(x)
        error('Simulator returned empty x.');
    end

    if isempty(u)
        warning('Simulator returned empty u.');
    end

    % --- normalize shapes ---
    t = t(:);

    % x should be N x (2*nq)
    if size(x,1) == 2*nq && size(x,2) ~= 2*nq
        x = x.';
    end

    % u should be N x nu
    if ~isempty(u)
        if size(u,1) == nu && size(u,2) ~= nu
            u = u.';
        end
    end

    % --- validate shapes ---
    if isempty(x) || size(x,2) ~= 2*nq
        error('Simulator returned x with invalid size [%d %d], expected N x %d.', ...
            size(x,1), size(x,2), 2*nq);
    end

    if isempty(u)
        error('Simulator returned empty u.');
    end

    if size(u,2) ~= nu
        error('Simulator returned u with invalid size [%d %d], expected N x %d.', ...
            size(u,1), size(u,2), nu);
    end

    if size(u,1) ~= size(x,1)
        error('x and u have inconsistent lengths: size(x,1)=%d, size(u,1)=%d.', ...
            size(x,1), size(u,1));
    end

    if any(isnan(x(:))) || any(isnan(u(:))) || length(t) < 3
        cineq = 1e6*ones(200,1);
        ceq   = 1e6*ones(30,1);
        return;
    end

    qStart  = x(1,   1:nq).';
    dqStart = x(1,   nq+1:2*nq).';
    qEnd    = x(end, 1:nq).';
    dqEnd   = x(end, nq+1:2*nq).';

    %% 1. Impact + relabelling
%     xPlus      = rabbit_impact_map([qEnd; dqEnd], p);
    xRelabeled = rabbit_reset_map([qEnd; dqEnd], params);
    qRel  = xRelabeled(1:nq);
    dqRel = xRelabeled(nq+1:2*nq);

    %% 2. Periodicity
    ceq = [ceq; qRel - qStart; dqRel - dqStart];

    %% 3. Speed
    [~,swing_foot_s,~,~,~,~]    = rabbit_kinematics(qStart, p);
    [~,swing_foot_e,~,~,~,~]    = rabbit_kinematics(qEnd,   p);
    stepLen = swing_foot_e(1) - swing_foot_s(1);
    ceq = [ceq; stepLen/T - opt.v_des];

    %% 4. Foot height at impact
    ceq = [ceq; swing_foot_e(2)];

    %% 5. HZD invariance
    [y0,    dy0   ] = hzd_virtualConstraints(qStart, dqStart, CP, model, opt);
    [yPlus, dyPlus] = hzd_virtualConstraints(qRel,   dqRel,  CP, model, opt);
    ceq = [ceq; y0; dy0; yPlus; dyPlus];

    %% 6. Path inequality constraints
    N = size(x,1);
    footClr   = zeros(N,1);
    hipViol   = zeros(N,1);
    torqueHi  = zeros(N, nu);
    torqueLo  = zeros(N, nu);
    jointHi   = zeros(N, nq);
    jointLo   = zeros(N, nq);
    kneeViol  = zeros(N,2);

    for k = 1:N
        qk = x(k,1:nq).';
        uk = u(k,:).';

        kin = rabbit_kinematics(qk, params);

        footClr(k)     = -kin.swingFoot(2);
        hipViol(k)     = opt.hipHeightMin - kin.hip(2);
        jointHi(k,:)   = qk.' - opt.qMax.';
        jointLo(k,:)   = opt.qMin.' - qk.';
        kneeViol(k,1)  = qk(5);
        kneeViol(k,2)  = qk(7);
        torqueHi(k,:)  = uk.' - opt.uMax;
        torqueLo(k,:)  = opt.uMin - uk.';
    end

    cineq = [footClr; hipViol; jointHi(:); jointLo(:); ...
        kneeViol(:); torqueHi(:); torqueLo(:)];

catch ME
    fprintf(2, '\n[hzd_constraintsHZD ERROR]\n');
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));

    keyboard

    c = 1e6;
    ceq = 1e6;
end
end
