function Zp = ch3_zd_point(theta, alpha, p)
%CH3_ZD_POINT  The zero dynamics surface Z at one value of the phase variable.
%
%   Zp = ch3_zd_point(theta, alpha, p)
%
% Z = {x : y = 0, ydot = 0} intersected with the stance contact constraint is a
% TWO-dimensional manifold, coordinatized by (theta, thetadot).  This function
% returns the point of Z at a given theta, together with the coefficients of the
% scalar second-order ODE the dynamics restrict to there.
%
% CONFIGURATION IS CLOSED FORM, NOT A NEWTON SOLVE.  The three defining
% conditions are triangular in this coordinate ordering:
%
%   y = 0             ->  q(iact)  = yd(s(theta), alpha)     (4 joints)
%   c_theta q = theta ->  q(3)     = (theta - c(4:7) q(4:7)) / c(3)
%   P_st(q) = 0       ->  q(1:2)   = -(dP/d[px pz]) \ P_st(q)|_{px=pz=0}
%
% The last step is one linear solve, and it is EXACT rather than a first Newton
% step, because the stance-foot position is affine in the floating base: moving
% (px,pz) translates the whole robot rigidly, so d P_st / d[px pz] is a constant
% 2x2 block of J_st.  ch3_test_hzd asserts the residual is at machine epsilon,
% which is the cheap way to keep that claim honest.
%
% VELOCITY IS ONE LINEAR SOLVE.  On Z, ydot = Jy qdot = 0 and the stance foot is
% fixed, so qdot is pinned up to the single free rate thetadot:
%
%       [ Jy    ]          [ 0 ]
%       [ c_theta ] qdot = [ 1 ] thetadot   ->   qdot = v_z(theta) thetadot
%       [ J_st  ]          [ 0 ]
%
% THE RESTRICTED DYNAMICS.  Substituting q = q_z(theta), qdot = v_z thetadot,
% and therefore qddot = v_z thetaddot + v_z' thetadot^2, into the constrained
% equations of motion and projecting onto the direction that sees NEITHER the
% actuators NOR the contact wrench gives a scalar ODE
%
%       a(theta) thetaddot + b(theta) thetadot^2 + c(theta) = 0.
%
% The projection row w is the one satisfying w B = 0 and w J_st' = 0.  Since
% B = [0_{3x4}; I_4], w B = 0 forces w(4:7) = 0, and what remains is the null
% vector of the 2x3 block J_st(:,1:3) -- computed here as the cross product of
% its two rows, which for a 2x3 matrix IS the null direction, with no svd and no
% tolerance to choose.
%
% WHY THE SCALE OF w DOES NOT MATTER.  Rescaling w by any function of theta
% scales a, b and c by the SAME factor, and everything downstream
% (ch3_zero_dynamics) uses only the ratios b/a and c/a.  So the normalization is
% free; the sign is fixed here by requiring a > 0, which makes a the effective
% inertia about the contact point and keeps the sign of sigma continuous.
%
% Inputs
%   theta : scalar phase variable
%   alpha : ny x n_ctrl Bezier coefficients
%   p     : parameter struct
%
% Outputs
%   Zp : struct with
%          .q .v      7x1 configuration on Z, and dq/dtheta (so qdot = v*thetadot)
%          .a .b .c   scalar coefficients of the restricted ODE
%          .w         1x7 projection row actually used
%          .s         the phase s(theta)
%          .res_foot  |P_st(q)|, the closed-form solve's residual
%
% See also CH3_ZERO_DYNAMICS, CH3_OUTPUTS, CH3_CONTROL_AFFINE.

nq = p.nq;
c_row = p.c_theta(:).';
dth   = p.theta_plus - p.theta_minus;

% --- configuration on Z ---------------------------------------------------
s = (theta - p.theta_minus) / dth;
[yd, dyd] = ch3_yd(alpha, s, p);

q = zeros(nq, 1);
q(p.iact) = yd;

% c_theta q = theta, solved for the one coordinate c_theta weights outside the
% actuated set (the torso pitch).  Written generically rather than as q(3) so
% that changing p.c_theta or p.iact cannot silently desynchronize this file.
i_free = find(c_row ~= 0 & ~ismember(1:nq, p.iact(:).'));
if numel(i_free) ~= 1
    error('ch3_zd_point:phaseStructure', ...
          ['c_theta must weight exactly one unactuated coordinate for the ' ...
           'configuration on Z to be closed form (found %d).'], numel(i_free));
end
q(i_free) = (theta - c_row(p.iact) * q(p.iact)) / c_row(i_free);

% floating base from the stance-foot gauge P_st(q) = 0.  Affine, so exact.
J_base = J_st(q);
J_base = J_base(:, 1:2);
q(1:2) = -(J_base \ P_st(q));

res_foot = norm(P_st(q), inf);

% --- velocity direction on Z ----------------------------------------------
ds_dq = c_row / dth;
Jy    = p.H - dyd * ds_dq;

Amat = [Jy; c_row; J_st(q)];
evec = [zeros(p.ny,1); 1; zeros(2,1)];

if rcond(Amat) < 1e-12
    error('ch3_zd_point:degenerate', ...
          ['The map from thetadot to qdot on Z is singular at theta = %.4f. ' ...
           'Z is not a graph over theta here.'], theta);
end
v = Amat \ evec;

% --- projection onto the unactuated, unconstrained direction --------------
Jst = J_st(q);
w3  = cross(Jst(1,1:3), Jst(2,1:3));
w   = [w3(:).', zeros(1, nq-3)];

Mq = M(q);
Vq = V([q; v]);            % Coriolis vector at thetadot = 1; quadratic in qdot
Gq = G(q);

a = w * Mq * v;

% v' = dv_z/dtheta by central difference.  v_z is smooth in theta (it is the
% solution of a linear system with smoothly varying entries), so a central
% difference is second-order accurate and needs no derivative of ch3_yd beyond
% the analytic dyd that ch3_bezier already provides.
hstep = 1e-6 * max(1, abs(dth));
vp = (vel_only(theta + hstep, alpha, p) - vel_only(theta - hstep, alpha, p)) / (2*hstep);

b = w * Mq * vp + w * Vq;
cc = w * Gq;

if a < 0
    a = -a;  b = -b;  cc = -cc;  w = -w;      % see the scale note in the header
end

if abs(a) < 1e-12
    error('ch3_zd_point:noInertia', ...
          ['The effective inertia a(theta) vanishes at theta = %.4f, so the ' ...
           'zero dynamics is not well posed there.'], theta);
end

Zp = struct('q', q, 'v', v, 'a', a, 'b', b, 'c', cc, ...
            'w', w, 's', s, 'res_foot', res_foot);

end

% ---------------------------------------------------------------------------
function v = vel_only(theta, alpha, p)
%VEL_ONLY  The v_z(theta) half of the computation, for the finite difference.
nq = p.nq;
c_row = p.c_theta(:).';
dth   = p.theta_plus - p.theta_minus;

s = (theta - p.theta_minus) / dth;
[yd, dyd] = ch3_yd(alpha, s, p);

q = zeros(nq, 1);
q(p.iact) = yd;
i_free = find(c_row ~= 0 & ~ismember(1:nq, p.iact(:).'));
q(i_free) = (theta - c_row(p.iact) * q(p.iact)) / c_row(i_free);
J_base = J_st(q);
q(1:2) = -(J_base(:,1:2) \ P_st(q));

Jy = p.H - dyd * (c_row / dth);
v  = [Jy; c_row; J_st(q)] \ [zeros(p.ny,1); 1; zeros(2,1)];

end
