%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRADIENT OF THE PHASE VARIABLE: d(theta)/dq = c, a CONSTANT row vector.
%
% Since theta_of_q is linear (theta = qt + q1 + 0.5*q2 = c*q), its Jacobian
% is just c -- exact, constant, and free of the singularities the previous
% atan2-based version had (it divided by |hip-foot|^2 and needed J_st).
%
% Keep this in sync with theta_of_q.m: they must be the same c.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = dtheta_dq_of(~)
    g = [0, 0, 1, 1, 0.5, 0, 0];
end
