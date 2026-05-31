function u = rabbit_controller_bspline(t, x, params)
% PD controller for RABBIT tracking B-spline virtual constraints
%
% Phase variable: torso angle qt (unactuated)
% s = (qt - theta0) / (thetaf - theta0)
%
% Virtual constraints: y = qa - hd(s)
%
% Inputs:
%   t:      current time
%   x:      state [q; dq] = [x;z;qt;q1;q2;q3;q4; dx;dz;dqt;dq1;dq2;dq3;dq4]
%   params: structure containing:
%       - ControlPoints: (n+1) x 4 B-spline control points
%       - n:             number of data points
%       - p:             B-spline degree
%       - theta0:        torso angle at start of step
%       - thetaf:        torso angle at end of step
%       - Kp, Kd:        PD gain matrices

    nq = 7;
    
    % Extract state
    q  = x(1:nq);
    dq = x(nq+1:end);
    
    % Extract parameters
    ControlPoints = params.ControlPoints;
    n     = params.n;
    p     = params.p;
    theta0 = params.theta0;
    thetaf = params.thetaf;
    Kp    = params.Kp;
    Kd    = params.Kd;
    
    % Compute phase variable from torso angle
    qt = q(3);
    s = (qt - theta0) / (thetaf - theta0);
    s = min(max(s, 0), 1);  % clamp to [0,1]
    
    % Compute phase derivative
    dqt = dq(3);
    ds_dt = dqt / (thetaf - theta0);
    
    % Evaluate desired trajectory and its derivative
    [hd, dhd_ds] = desired_gait_bspline(s, ControlPoints, n, p);
    
    % Virtual constraints and their derivatives
    qa  = q(4:7);             % actuated joints [q1;q2;q3;q4]
    dqa = dq(4:7);
    
    y  = qa  - hd;                  % position error
    dy = dqa - dhd_ds * ds_dt;      % velocity error
    
    % PD control law
    u = -Kp * y - Kd * dy;
    
end