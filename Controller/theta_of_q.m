%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE VARIABLE: absolute stance leg angle, LINEAR in q.
%
%   theta = qt + q1 + 0.5*q2  =  c*q,   c = [0 0 1 1 0.5 0 0]
%
% This is Hypothesis HH6 of Westervelt et al., "Feedback Control of Dynamic
% Bipedal Robot Locomotion" (Sec. 6.4.1): the phase variable must be an
% ABSOLUTE ANGLE that is LINEAR in the configuration, theta = c*q. The whole
% HZD construction relies on it (theta cyclic => the inertia matrix reduces,
% eqs. 6.55-6.68), and it is what makes the decoupling-matrix analysis valid.
%
% This is NOT an approximation of the previous atan2 formulation -- it is
% EXACTLY the same function with the branch cut removed. Because both leg
% links have length 1/2, the stance hip->foot vector satisfies
%
%   hip - foot = cos(q2/2) * [ sin(qt+q1+q2/2) ; cos(qt+q1+q2/2) ]
%
% so atan2(rel_x, rel_z) == qt + q1 + q2/2 whenever cos(q2/2) > 0, i.e. for
% every |q2| < pi (the knee never folds through 180 deg). Verified to 3e-13
% over thousands of random poses.
%
% Why it matters: atan2 WRAPS at +/-pi. Once the optimizer pushed the leg
% past that (theta ~ 2.6 rad) the phase variable jumped to -3.04 and the
% gait optimization chased a discontinuity. The linear form cannot wrap,
% and its gradient is the constant c (see dtheta_dq_of.m) instead of a
% nonlinear expression that divides by |hip-foot|^2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function th = theta_of_q(q)
    th = q(3) + q(4) + 0.5*q(5);
end
