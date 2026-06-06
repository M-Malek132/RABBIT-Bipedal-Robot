function h0 = hzd_h0_outputs(q, model)
%HZD_H0_OUTPUTS  Select the ny controlled outputs from q.
%
%  For RABBIT with q = [px, pz, qt, q1, q2, q3, q4]':
%    - px, pz are unactuated (floating base)
%    - qt is unactuated (torso)
%    - q1, q2, q3, q4 are the 4 actuated joint angles
%
%  Therefore h0 = [q1; q2; q3; q4] = q(4:7)

h0 = q(4:7);    % ny x 1  (ny = 4)
end
