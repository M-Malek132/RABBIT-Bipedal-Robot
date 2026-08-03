function pv = ch5_pend_pv(p)
%CH5_PEND_PV  Pack the pendulum parameters into the canonical 12-vector.
%
%   pv = ch5_pend_pv(p)
%
% ONE ORDERING, DECLARED ONCE. ch5_pendulum unpacks it, ch5_gen_pendulum
% differentiates against it, and the generated code in Model/Generated takes it
% verbatim. Reordering this vector silently invalidates every committed
% generated file, so if it ever changes, regenerate -- ch5_test_model compares
% the generated derivatives against the live model and will catch a mismatch,
% but only if the test is run.
%
%   1  m1     2  m2      link masses               [kg]
%   3  l1     4  l2      link lengths              [m]
%   5  lc1    6  lc2     COM distance along link   [m]
%   7  I1     8  I2      link inertia about COM    [kg m^2]
%   9  g                 gravity                   [m/s^2]
%  10  Jm                motor inertia             [kg m^2]
%  11  k                 motor/joint stiffness     [Nm/rad]
%  12  xi                joint damping             [Nm s/rad]
%
% See also CH5_PENDULUM, CH5_GEN_PENDULUM.

q = p.plant;
pv = [q.m(1); q.m(2); q.l(1); q.l(2); q.lc(1); q.lc(2); ...
      q.I(1);  q.I(2); q.g;   q.Jm;   q.k;     q.xi];

end
