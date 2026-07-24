function [x_plus, impulse] = rabbit_impact_map(x_minus)
%RABBIT_IMPACT_MAP  Rigid, instantaneous swing-foot impact (velocity reset).
%
%   [x_plus, impulse] = rabbit_impact_map(x_minus)
%
% Applies the plastic-impact model of Westervelt et al. (Sec. 3.4.2 / eq. 3.20):
% when the swing foot strikes the ground it sticks (no rebound, no slip), which
% imposes the holonomic constraint J_sw*dq_plus = 0. Solving the momentum
% balance for the impulsive contact force gives the KKT system
%
%       [ M  -J_sw' ] [ dq_plus ]   [ M*dq_minus ]
%       [ J_sw   0   ] [ lambda  ] = [     0      ]
%
% The configuration q is CONTINUOUS through impact (only velocities jump); the
% impulsive force lives entirely in the velocity reset.
%
% IMPORTANT — this is the PRE-relabel map. J_sw and M are evaluated at the
% PRE-IMPACT configuration q = x_minus(1:7); the stance/swing leg RELABELING
% happens afterward in rabbit_reset_map, which every caller applies as
% rabbit_reset_map(rabbit_impact_map(x_minus)). Do not swap the order.
%
% Inputs
%   x_minus : 14x1 pre-impact state [q; dq] (q in generalized coords, rad/m)
% Outputs
%   x_plus  : 14x1 post-impact state [q; dq_plus]  (q unchanged)
%   impulse : 2x1 contact impulse [Ix; Iz] at the swing foot [Ns], world frame
%             (up-positive z). Optional 2nd output; PD callers use only x_plus.

    q  = x_minus(1:7);
    dq = x_minus(8:14);

    D = M(q);
    J = J_sw(q);   % swing-foot Jacobian at the PRE-impact configuration

    % System of equations:
    % [ M  -J' ] [ dq_plus ] = [ M * dq_minus ]
    % [ J   0  ] [ lambda  ]   [ 0            ]

    A = [D, -J';
         J,  zeros(size(J,1))];
    b = [D*dq;
         zeros(size(J,1),1)];

    sol = A\b;
    dq_plus = sol(1:7);

    % lambda here is the CONTACT IMPULSE [Ns] at the swing foot (world [x;z]).
    % Returned as an optional 2nd output; existing callers use only x_plus.
    impulse = sol(8:end);

    % Note: q does not change during impact, only dq
    x_plus = [q; dq_plus];
end
