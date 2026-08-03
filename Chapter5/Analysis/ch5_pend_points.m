function P = ch5_pend_points(theta, p)
%CH5_PEND_POINTS  Joint positions of the 2-link pendulum, world frame.
%
%   P = ch5_pend_points(theta, p)
%
% theta is 2x1 (or 2xN). Returns a 2 x 3 x N array of [base, elbow, tip], with
% x to the right and y UP-POSITIVE -- the same convention as the rest of the
% repository's kinematic helpers (see the note in README about world-frame Z).
%
% theta1 is measured from the DOWNWARD vertical and theta2 relative to link 1,
% so theta = 0 hangs straight down and theta1 = +-pi points straight up. The
% tip's y coordinate is the py2 that ch5_barrier constrains.
%
% See also CH5_PENDULUM, CH5_ANIMATE_PENDULUM, CH5_PLOT_PENDULUM.

l1 = p.plant.l(1);
l2 = p.plant.l(2);

th1 = theta(1,:);
th2 = theta(2,:);
N   = numel(th1);

base  = zeros(2, N);
elbow = [ l1*sin(th1); -l1*cos(th1) ];
tip   = elbow + [ l2*sin(th1+th2); -l2*cos(th1+th2) ];

P = cat(3, base, elbow, tip);      % 2 x N x 3
P = permute(P, [1 3 2]);           % 2 x 3 x N

end
