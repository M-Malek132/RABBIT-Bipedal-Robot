function x0 = ch5_x0(p)
%CH5_X0  The initial state of the selected study.
%
%   x0 = ch5_x0(p)
%
% springmass  everything at the origin, at rest. The target x3d = 3 m is then
%             three metres away through six integrators, which is what makes
%             Fig. 5.3 a real transient rather than a regulation problem.
%
% pendulum    theta = [-pi; 0] with the motors UNLOADED, thetam = theta, and
%             everything at rest.
%
%             That is a genuine equilibrium, not merely a starting pose: with
%             theta1 = -pi and theta2 = 0 the arm is inverted and colinear, so
%             gravity produces no torque; with thetam = theta the springs are
%             unstretched, so (5.37) produces none either. xdot(0) = 0 exactly,
%             which ch5_test_model asserts. The maneuver therefore starts from
%             rest at an UNSTABLE equilibrium -- the controller has to do all
%             the work, and none of the motion is inherited from a badly posed
%             initial condition.
%
% See also CH5_SYSTEM, CH5_SIMULATE.

switch lower(p.system)

    case 'springmass'
        x0 = p.plant.x0(:);

    case 'pendulum'
        th0 = p.plant.theta0(:);
        x0  = [th0; th0; zeros(4,1)];    % [theta; thetam; thetadot; thetamdot]

    otherwise
        error('ch5_x0:unknownSystem', 'Unknown p.system "%s".', p.system);
end

end
