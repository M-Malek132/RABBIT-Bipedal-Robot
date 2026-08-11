function k = ch6_kin(x, aux, p, point)
%CH6_KIN  A tracked point, relative to the stance foot, to SECOND order in u.
%
%   k = ch6_kin(x, aux, p, point)          point = 'swing' | 'head'
%   k = ch6_kin(x, [],  p, point)          computes aux itself
%
% Every barrier in Chapter 6 is a function of ONE two-dimensional point
% measured from the stance foot:
%
%   'swing'   r = (l_f, h_f)   swing foot     Sections 6.1.3, 6.2, 6.4
%   'head'    r = (l_H, h_H)   torso top      Section 6.1.2
%
% and every barrier row needs that point differentiated TWICE along
% xdot = f + g u, because a position constraint has relative degree 2. This
% function returns exactly the pieces that makes possible and nothing else:
%
%       rddot = Lf2r + LgLfr * u
%
% -------------------------------------------------- why relative, not absolute
% l_f and h_f are DIFFERENCES, P_sw - P_st, not P_sw alone. Three reasons, and
% only the first is obvious:
%
%   1. It is what Fig. 6.3 draws. The circles O1, O2 are centred on the stance
%      foot, so the constraint is naturally written in the stance frame.
%   2. The stance foot DRIFTS. The contact constraint is imposed at the
%      acceleration level, so ode45 lets the pinned foot wander slowly within a
%      step (the repo README calls this out; ch3_report measures ~9e-07 m on the
%      reference gait). An absolute h_f would inherit that drift as a bias in
%      the barrier; a relative one cancels it exactly.
%   3. It survives the re-plant in Delta. ch3_impact shifts pz so the NEW stance
%      foot sits at z = 0, which moves every absolute height by a small amount
%      at every impact. Differences do not move.
%
% ---------------------------------------------------- why no kinematic Hessian
% The barriers are nonlinear in r -- circle distances, mostly -- so gddot needs
% a second derivative somewhere. It does NOT need the second derivative of the
% KINEMATICS. Writing g = G(r) and differentiating twice,
%
%       gddot = dG/dr * rddot  +  rdot' * d2G/dr2 * rdot
%
% puts all the curvature in d2G/dr2, a 2x2 matrix of the BARRIER, which is
% closed form for every barrier in the chapter. The kinematics only ever appear
% through rddot, which the generated Jdotdq_* files already give exactly. So
% there is no symbolic Hessian of P_sw anywhere in Chapter 6, and nothing here
% is finite-differenced.
%
% Inputs
%   x     : 14x1 state
%   aux   : the third output of ch3_control_affine, or [] to compute it
%   p     : parameter struct (uses p.nq, p.nu)
%   point : 'swing' or 'head'
%
% Output
%   k : struct with
%         .r      2x1   [horizontal; vertical], relative to the stance foot [m]
%         .rdot   2x1   drdot/dt
%         .Jr     2xnq  rdot = Jr dq
%         .Jdr    2x1   Jrdot dq
%         .Lf2r   2x1   Jr*ddq_drift + Jdr        (the u-free part of rddot)
%         .LgLfr  2xnu  Jr*ddq_in                 (the coefficient of u)
%         .aux    the aux struct, so callers reuse one KKT solve
%
% See also CH6_BARRIER, CH3_CONTROL_AFFINE, CH6_BAR_STONES.

nq = p.nq;
q  = x(1:nq);
dq = x(nq+1:2*nq);

if isempty(aux)
    [~, ~, aux] = ch3_control_affine(x, p);
end

% --- stance foot: the origin of every Chapter-6 coordinate -----------------
J_s  = J_st(q);              % 2 x nq
Jd_s = Jdotdq_st(q, dq);     % 2 x 1
r_s  = P_st(q);

switch lower(point)

    case 'swing'
        r_p  = P_sw(q);
        J_p  = J_sw(q);
        Jd_p = Jdotdq_sw(q, dq);

    case 'head'
        % Torso top, from Tt(q)*[0; -L_T; 0; 1]:
        %
        %       x_H = px + L_T sin(qt)
        %       z_H = L_T cos(qt) - pz
        %
        % Hand-written rather than generated because it is three lines and its
        % Jacobian is exact by inspection -- the only q-dependence is through
        % qt, so J and Jdot*dq below are complete, not truncated.
        %
        % L_T must match ch3_body_points and the torso length baked into
        % Dynamics/rabbit_energy_model_generalized_Lagrange.m, or the barrier
        % protects a head the simulation does not have.
        L_T = 0.75;
        qt  = q(3);
        dqt = dq(3);

        r_p = [q(1) + L_T*sin(qt); ...
               L_T*cos(qt) - q(2)];

        J_p = zeros(2, nq);
        J_p(1, 1) =  1;   J_p(1, 3) =  L_T*cos(qt);
        J_p(2, 2) = -1;   J_p(2, 3) = -L_T*sin(qt);

        Jd_p = [-L_T*sin(qt) * dqt^2; ...
                -L_T*cos(qt) * dqt^2];

    otherwise
        error('ch6_kin:point', ...
              'Unknown tracked point "%s" (expected swing|head).', point);
end

Jr  = J_p  - J_s;
Jdr = Jd_p - Jd_s;

k = struct();
k.r     = r_p - r_s;
k.Jr    = Jr;
k.Jdr   = Jdr;
k.rdot  = Jr * dq;
k.Lf2r  = Jr * aux.ddq_drift + Jdr;
k.LgLfr = Jr * aux.ddq_in;
k.aux   = aux;

end
