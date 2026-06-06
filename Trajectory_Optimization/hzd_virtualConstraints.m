function [y, dy, Jy] = hzd_virtualConstraints(q, dq, CP, model, opt)
%HZD_VIRTUALCONSTRAINTS  Compute virtual constraint outputs using B-spline.
%
%  y  = h0(q) - hd(s(q))                 (4×1)
%  dy = Jy * dq                           (4×1)
%  Jy = dh0/dq - dhd/dtheta * dtheta/dq  (4×7)
%
%  h0(q)    = q(4:7)        actual actuated joint angles
%  hd(s)    = CP' * N(s)    desired angles from B-spline
%  theta    = q(1) = px     phase variable (monotone)
%  s        = (theta - thetaStart) / (thetaEnd - thetaStart)
%
%  Consistent with BSplineTrajectory.virtual_constraint, but also
%  returns the full output Jacobian Jy needed by the HZD controller.

% --- Actual outputs and their Jacobian ---
h0  = q(4:7);                         % 4×1
Jh0 = hzd_jacobian_h0(q, model);      % 4×7  = [0_{4x3}, I_{4x4}]

% --- Phase variable ---
[theta, dtheta_dq] = hzd_phaseVariable(q, model);   % scalar, 1×7

% --- Desired outputs from B-spline ---
[hd, dhd_dtheta] = hzd_evalBSpline(theta, CP, opt);  % both 4×1

% --- Virtual constraint error ---
y = h0 - hd;                           % 4×1

% --- Output Jacobian: dy/dq ---
% Jy = Jh0  -  dhd_dtheta * dtheta_dq
%     (4×7)    (4×1)        (1×7)
Jy = Jh0 - dhd_dtheta * dtheta_dq;    % 4×7

% --- Time derivative of y ---
dy = Jy * dq;                          % 4×1
end
