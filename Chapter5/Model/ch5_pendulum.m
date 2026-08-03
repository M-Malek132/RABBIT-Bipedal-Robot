function d = ch5_pendulum(x, pv)
%CH5_PENDULUM  The 2-link pendulum with elastic actuators, Fig. 5.2, (5.36)-(5.37).
%
%   d = ch5_pendulum(x, pv)
%
% SYMBOLIC-SAFE BY CONSTRUCTION. Every operation below is plain algebra -- no
% branching on values, no backslash, no indexing that depends on data -- so
% this same function is what ch5_gen_pendulum differentiates with `syms` in
% place of numbers. There is therefore exactly ONE statement of these dynamics
% in the repository, and the generated Lie derivatives cannot drift away from
% the model the simulation integrates. That is worth more than the small cost
% of writing D^-1 out longhand instead of calling backslash.
%
% ------------------------------------------------------------------ the model
% Coordinates  q = [theta1; theta2; thetam1; thetam2],  x = [q; qdot].
% theta1 from the DOWNWARD vertical, theta2 relative to link 1, so theta = 0 is
% the arm hanging straight down and theta1 = +-pi is straight up.
%
% Link side, a standard planar 2R arm driven by the JOINT torque of (5.37):
%
%       D(theta) thetaddot + C(theta,thetadot) thetadot + Gv(theta) = u_joint
%       u_joint = -k(theta - thetam) - xi thetadot                      (5.37)
%
% Motor side, (5.36):
%
%       Jm thetamddot = k(theta - thetam) + tau
%
% ------------------------------------------------------- where the 4 comes from
% tau appears in thetamddot. It reaches thetaddot only through u_joint, which
% depends on thetam, NOT on thetamddot. So differentiate the output theta:
%
%       theta^(2) = D^-1(...)                    thetam appears, tau does not
%       theta^(3)                                thetamdot appears, tau does not
%       theta^(4)                                thetamddot appears -> TAU
%
% RELATIVE DEGREE 4, and with two outputs on eight states the arm is again
% fully linearized with no zero dynamics. Note what the elasticity did: a rigid
% 2R arm driven directly by tau is relative degree 2, and the series spring is
% what pushes it to 4. This is the chapter's stated motivation -- "mechanical
% systems with springs or series elastic actuators that are often seen in many
% dynamic walking robots" -- rather than a contrived example.
%
% Note also that g(x) is CONSTANT here: tau enters only the motor rows, through
% 1/Jm. The decoupling matrix L_g L_f^3 y = D^-1 k / Jm is not constant, but it
% is singular only where D is, i.e. never.
%
% Inputs
%   x  : 8x1 state (may be symbolic)
%   pv : 12x1 packed parameters, see ch5_pend_pv
%
% Output
%   d : struct with
%         .theta .thetam .dtheta .dthetam
%         .D .Cm .Gv .u_joint      the (5.37) joint torque actually applied
%         .thetaddot .thetamddot_drift
%         .f  8x1  drift              xdot = f(x) + g*tau
%         .g  8x2  input map (constant)
%         .py .Jpy                    end-effector height and its 1x2 gradient
%
% See also CH5_PEND_PV, CH5_GEN_PENDULUM, CH5_CONTROL_AFFINE.

m1 = pv(1);  m2 = pv(2);
l1 = pv(3);  l2 = pv(4);
lc1= pv(5);  lc2= pv(6);
I1 = pv(7);  I2 = pv(8);
gr = pv(9);
Jm = pv(10); k = pv(11); xi = pv(12);

th1  = x(1);  th2  = x(2);
thm1 = x(3);  thm2 = x(4);
dth1 = x(5);  dth2 = x(6);
dthm1= x(7);  dthm2= x(8);

theta   = [th1;  th2];
thetam  = [thm1; thm2];
dtheta  = [dth1; dth2];
dthetam = [dthm1; dthm2];

s2 = sin(th2);
c2 = cos(th2);

%% ------------------------------------------------------------ mass matrix
D11 = I1 + I2 + m1*lc1^2 + m2*(l1^2 + lc2^2 + 2*l1*lc2*c2);
D12 = I2 + m2*(lc2^2 + l1*lc2*c2);
D22 = I2 + m2*lc2^2;
D   = [D11, D12; D12, D22];

%% ---------------------------------------------------- Coriolis / centrifugal
hc  = -m2*l1*lc2*s2;
Cm  = [hc*dth2, hc*(dth1 + dth2); ...
       -hc*dth1, 0*hc];               % 0*hc keeps the class symbolic-safe

%% -------------------------------------------------------------- gravity
% From U = -m1 g lc1 cos(th1) - m2 g (l1 cos(th1) + lc2 cos(th1+th2)), so Gv
% vanishes at theta = 0 (hanging) and at th1 = pi, th2 = 0 (inverted) -- the
% latter is why x0 of Fig. 5.4 is a genuine equilibrium and the maneuver starts
% from rest rather than from a transient.
s12 = sin(th1 + th2);
G1  = (m1*lc1 + m2*l1)*gr*sin(th1) + m2*lc2*gr*s12;
G2  = m2*lc2*gr*s12;
Gv  = [G1; G2];

%% ------------------------------------------------- the two coupled subsystems
u_joint    = -k*(theta - thetam) - xi*dtheta;              % (5.37)
rhs        = u_joint - Cm*dtheta - Gv;

% Explicit 2x2 inverse rather than backslash: exact, fast, and the only form
% that survives symbolic differentiation without expression blowup.
detD       = D11*D22 - D12*D12;
Dinv       = [D22, -D12; -D12, D11] / detD;
thetaddot  = Dinv * rhs;

thetamddot_drift = (k/Jm) * (theta - thetam);              % (5.36) without tau

%% ------------------------------------------------------- control-affine form
f = [dtheta; dthetam; thetaddot; thetamddot_drift];
g = [zeros(6,2); eye(2)/Jm];

%% ------------------------------------------------------------ end effector
% py2 of Fig. 5.2, world frame, up-positive. Depends on theta only -- which is
% what gives the constraint the same relative degree as the outputs.
c1  = cos(th1);
c12 = cos(th1 + th2);
py  = -l1*c1 - l2*c12;

% 1x2 gradient dpy/dtheta. Worth naming: because the outputs ARE theta, this
% row turns out to be exactly the coefficient of mu in the VIOL barrier row
% (see ch5_ctrl_ecbf_clf_qp), so it is the thing that decides whether the
% safety row is controllable at a given configuration. It vanishes iff
% sin(th1+th2) = 0 and sin(th1) = 0, i.e. the arm fully extended vertically.
s1   = sin(th1);
Jpy  = [l1*s1 + l2*s12, l2*s12];

d = struct('theta', theta, 'thetam', thetam, ...
           'dtheta', dtheta, 'dthetam', dthetam, ...
           'D', D, 'Cm', Cm, 'Gv', Gv, 'u_joint', u_joint, ...
           'thetaddot', thetaddot, 'thetamddot_drift', thetamddot_drift, ...
           'f', f, 'g', g, 'py', py, 'Jpy', Jpy);

end
