function [theta, dtheta_dq] = hzd_phaseVariable(q, model)
%HZD_PHASEVARIABLE Posture-based phase variable for RABBIT.
%
% Coordinates:
%   q = [px, pz, qt, q1, q2, q3, q4]'
%
% Phase:
%   theta = -qt - q1 - 0.5*q2
%
% MATLAB indexing:
%   qt = q(3), q1 = q(4), q2 = q(5)

    nq = model.nq;
    % assumption of equal length for shank and thigh
    theta = +q(3) + q(4) + 0.5*q(5);

    dtheta_dq = zeros(1, nq);
    dtheta_dq(3) = +1.0;
    dtheta_dq(4) = +1.0;
    dtheta_dq(5) = +0.5;
end
