function [Mq, Vv, Gv] = ch4_link_dynamics(q, dq, L, g0)
%CH4_LINK_DYNAMICS  M(q), V(q,dq), G(q) of the planar robot for ANY link parameters.
%
%   [Mq, Vv, Gv] = ch4_link_dynamics(q, dq, L)
%   [Mq, Vv, Gv] = ch4_link_dynamics(q, dq, L, g0)
%
% The generated M.m / V.m / G.m have today's masses and inertias baked in, so
% they can express a uniformly scaled robot (every term times s) and a point
% load at the hip, and nothing else. A STRUCTURED uncertainty set -- each link
% heavier or lighter on its own, its COM shifted, its inertia off -- needs the
% equations of motion as a function of the link parameters. This is that
% function, numerically and exactly (no finite differences), for the same
% kinematic chain the Lagrangian generator uses:
%
%   hip          (x, -y)  in world (X, Z), Z up         [q = x y qt q1 q2 q3 q4]
%   torso        absolute angle phi = qt
%   stance leg   thigh qt + q1, shank qt + q1 + q2, knee at l_thigh
%   swing leg    thigh qt + q3, shank qt + q3 + q4
%
% A body at absolute angle phi has local axes u(phi) = (cos phi, -sin phi) and
% w(phi) = (-sin phi, -cos phi) in world (X, Z) -- the columns of the
% generator's transforms -- so a point at local (a, b) sits at a u + b w, with
% du/dphi = w and dw/dphi = -u. Each COM is the hip plus a sum of such levers,
% one per absolute angle it hangs from, which gives, per body i,
%
%   J_i      = dCOM_i/dq             (2 x 7)
%   abias_i  = -sum_k c_k phidot_k^2  (the COM acceleration at qddot = 0)
%   Jrot_i   = dphi_i/dq             (a constant row)
%
%   M = sum_i  m_i J_i' J_i + Izz_i Jrot_i' Jrot_i
%   V = sum_i  m_i J_i' abias_i          (planar bodies: no gyroscopic term)
%   G = sum_i  m_i g0 J_i(Z,:)'          (U = sum m g0 z_COM, z up)
%
% which is M qddot + V + G = B u, the generator's convention. With
% L = ch4_link_params() it reproduces M.m, V.m and G.m (ch4_test_model).
%
% Inputs
%   q, dq : 7x1 configuration and velocity (dq may be [] if V is not needed)
%   L     : 1x5 struct array as ch4_link_params returns, perturbed as desired
%   g0    : gravity [m/s^2], default 9.8062 (the generator's)
%
% See also CH4_LINK_PARAMS, CH4_UNCERTAINTY_SET, CH4_CONTROL_AFFINE.

if nargin < 4 || isempty(g0), g0 = 9.8062; end
q = q(:);
need_V = nargout > 1 && ~isempty(dq);
if need_V, dq = dq(:); end

% absolute angles as constant rows of q
Erow = zeros(5, 7);
Erow(1, 3)       = 1;          % torso
Erow(2, [3 4])   = 1;          % stance thigh
Erow(3, [3 4 5]) = 1;          % stance shank
Erow(4, [3 6])   = 1;          % swing thigh
Erow(5, [3 6 7]) = 1;          % swing shank
phi = Erow * q;
if need_V, phid = Erow * dq; end

Jhip = [1 0 0 0 0 0 0; ...
        0 -1 0 0 0 0 0];       % d(hip)/dq: hip = (x, -y)

% levers per body: rows [angle index, a, b]
lev = { [1, L(1).a, L(1).b]; ...
        [2, L(2).a, L(2).b]; ...
        [2, 0, L(2).l; 3, L(3).a, L(3).b]; ...
        [4, L(4).a, L(4).b]; ...
        [4, 0, L(4).l; 5, L(5).a, L(5).b] };

Mq = zeros(7);
Vv = zeros(7, 1);
Gv = zeros(7, 1);
for i = 1:5
    J  = Jhip;
    ab = [0; 0];
    T  = lev{i};
    for r = 1:size(T, 1)
        k = T(r, 1); A = T(r, 2); Bb = T(r, 3);
        f  = phi(k);
        uu = [cos(f); -sin(f)];
        ww = [-sin(f); -cos(f)];
        J  = J + (A * ww - Bb * uu) * Erow(k, :);
        if need_V
            ab = ab - (A * uu + Bb * ww) * phid(k)^2;
        end
    end
    Mq = Mq + L(i).m * (J.' * J) + L(i).J * (Erow(i, :).' * Erow(i, :));
    if need_V
        Vv = Vv + L(i).m * (J.' * ab);
    end
    Gv = Gv + L(i).m * g0 * J(2, :).';
end
Mq = (Mq + Mq.') / 2;
if ~need_V, Vv = []; end

end
