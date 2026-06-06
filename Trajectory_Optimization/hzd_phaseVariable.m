function [theta, dtheta_dq] = hzd_phaseVariable(q, model)
%HZD_PHASEVARIABLE  Phase variable for 7-DOF RABBIT — consistent with
%                   BSplineTrajectory.phase_variable.
%
%  q = [px, pz, qt, q1, q2, q3, q4]'   (7×1)
%
%  theta = q(1) = px  (forward base position)
%
%  This is monotonically increasing during a walking step and matches
%  exactly the phase variable used inside BSplineTrajectory:
%      theta = q(1);   ds = dq(1) / (thetaf - theta0)
%
%  opt.thetaStart and opt.thetaEnd are therefore px values [m].

theta       = q(1);
dtheta_dq   = [1, 0, 0, 0, 0, 0, 0];   % 1×7 row vector
end
