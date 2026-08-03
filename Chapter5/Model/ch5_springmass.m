function sm = ch5_springmass(p)
%CH5_SPRINGMASS  The serial spring-mass system of Fig. 5.1, eq (5.32)-(5.35).
%
%   sm = ch5_springmass(p)
%
% Three carts in a line, springs between neighbours, and the force acts ON THE
% FIRST CART ONLY:
%
%       m1 xddot1 = u + k(x2 - x1)
%       m2 xddot2 = k(x1 - x2) + k(x3 - x2)
%       m3 xddot3 = k(x2 - x3)
%
% With x = [x1 x2 x3 xdot1 xdot2 xdot3]' this is xdot = A x + B u, eq (5.35).
%
% ------------------------------------------- why this system is in Chapter 5
% Count the integrators between u and x3. The force has to travel cart 1 ->
% spring -> cart 2 -> spring -> cart 3, and each cart contributes two:
%
%       x3^(1)     xdot3                                        no u
%       x3^(2)     (k/m3)(x2 - x3)                              no u
%       x3^(3)     (k/m3)(xdot2 - xdot3)                        no u
%       x3^(4)     (k/m3)[(k/m2)(x1-2x2+x3) - (k/m3)(x2-x3)]    no u
%       x3^(5)     the same in velocities                       no u
%       x3^(6)     ... + (k/m3)(k/m2)(1/m1) u                   <- u at last
%
% RELATIVE DEGREE 6, on a 6-state system. So the output x3 linearizes the
% ENTIRE state -- there are no zero dynamics to argue about -- and the safety
% constraint x3 <= x3max, being a constraint on that same output, inherits
% relative degree 6 as well. A reciprocal CBF (Section 5.1) has literally
% nothing to say here: Lg(1/h) = 0 identically, so its barrier row does not
% contain the control at all. That is the chapter's motivating failure, made
% concrete, and ch5_test_qp asserts it.
%
% With the default unit parameters the relevant Lie derivatives come out as
%
%       L_g L_f^5 h  =  k^2/(m1 m2 m3)  =  1
%       L_f^6 y      =  -4 x1 + 9 x2 - 5 x3        (k = m_i = 1)
%
% and the second of those vanishes on x1 = x2 = x3, which is why the
% feedforward u_ff -> 0 at the target and the plotted force settles to zero
% rather than to a holding value.
%
% Output
%   sm : struct .A .B .n .nu  and .x0
%
% Cached: A and B depend only on (m, k), and this sits behind an ODE
% right-hand side.
%
% See also CH5_CONTROL_AFFINE, CH5_IO_LIN, CH5_BARRIER, CH5_SYSTEM.

persistent key val

m = p.plant.m(:).';
k = p.plant.k;

this = {m, k};
if ~isempty(key) && isequaln(key, this)
    sm = val;
    sm.x0 = p.plant.x0(:);
    return;
end

m1 = m(1); m2 = m(2); m3 = m(3);

% Stiffness of the free-free chain: row i is the force on cart i per unit
% displacement. Note the interior cart sees BOTH springs, hence the -2k.
Kc = [ -k,    k,    0 ; ...
        k,  -2*k,   k ; ...
        0,    k,   -k ];

Minv = diag(1 ./ [m1 m2 m3]);

A = [ zeros(3), eye(3) ; ...
      Minv*Kc,  zeros(3) ];

B = [ 0; 0; 0; 1/m1; 0; 0 ];

sm = struct('A', A, 'B', B, 'n', 6, 'nu', 1, 'x0', p.plant.x0(:));

key = this;
val = sm;

end
