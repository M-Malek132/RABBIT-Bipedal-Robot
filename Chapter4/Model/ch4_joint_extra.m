function tau = ch4_joint_extra(dq, S)
%CH4_JOINT_EXTRA  Joint torques the true actuators add to the command.
%
%   tau = ch4_joint_extra(dq, S)
%
% Friction and bias at the four actuated joints, in the joint (output) frame,
% so that the robot receives u + tau:
%
%   tau = -b_visc .* qdot - tau_coulomb .* tanh(qdot / coulomb_vel) + tau_bias
%
% with qdot = dq(4:7), the relative joint rates. The Coulomb term is smoothed
% by tanh so the plant stays Lipschitz for the integrator; at 0.05 rad/s it is
% within 1% of the sign function beyond 0.13 rad/s. S = [] gives zeros.
%
% See also CH4_STRUCTURED_PART, CH4_CONTROL_AFFINE.

if isempty(S)
    tau = zeros(4, 1);
    return;
end
w   = dq(4:7);
w   = w(:);
tau = -S.b_visc .* w - S.tau_coulomb .* tanh(w / S.coulomb_vel) + S.tau_bias;
end
