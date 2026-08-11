function io = ch5_io_lin(x, p)
%CH5_IO_LIN  Input-output linearization of the tracking outputs, eq (3.9)/(3.13).
%
%   io = ch5_io_lin(x, p)
%
% Differentiate the output r times until the control appears:
%
%       y^(r) = L_f^r y  +  L_g L_f^(r-1) y  u                      (5.21)
%
% and pre-invert the decoupling matrix,
%
%       u = u_ff + (L_g L_f^(r-1) y)^-1 mu,   u_ff = -(...)^-1 L_f^r y   (3.9)
%
% so that the transverse dynamics collapse to ny parallel r-integrator chains,
%
%       etadot = F eta + G mu,   eta = [y; ydot; ...; y^(r-1)].       (3.13)
%
% Both Chapter-5 plants have r*ny = nx, so eta is a change of coordinates on the
% WHOLE state and there are no zero dynamics. Nothing the controller does can
% be hidden in an unobserved internal mode -- a convenience of these two
% validation systems that Chapter 6's walking robot will not share.
%
% ------------------------------------------------------------ the springmass
% Linear, so every Lie derivative is a matrix power and is exact:
%
%       y      = Cy x - x3d,   Cy = e3'
%       y^(i)  = Cy A^i x           (a constant offset differentiates away)
%       L_f^6 y      = Cy A^6 x
%       L_g L_f^5 y  = Cy A^5 B = k^2/(m1 m2 m3)
%
% -------------------------------------------------------------- the pendulum
% Nonlinear, so the same quantities come from the committed symbolic
% derivation in Model/Generated. Because the outputs ARE theta, eta is just
% [theta - theta_d; thetadot; thetaddot; theta^(3)].
%
% Output
%   io : struct with
%          .y        ny x 1        the output itself
%          .eta      (r*ny) x 1    [y; ydot; ...; y^(r-1)]
%          .Lfry     ny x 1        L_f^r y
%          .LgLfr1y  ny x nu       L_g L_f^(r-1) y, the decoupling matrix
%          .Ainv     nu x ny       its inverse (or pinv if ill-conditioned)
%          .u_ff     nu x 1
%          .rcond    conditioning of the decoupling matrix
%
% See also CH5_CONTROL_AFFINE, CH5_BARRIER, CH5_CTRL_CLF_QP.

x = x(:);

switch lower(p.system)

    case 'springmass'
        sm = ch5_springmass(p);
        Cy = [0 0 1 0 0 0];                    % y tracks x3
        r  = p.sys.r;

        eta = zeros(r, 1);
        row = Cy;
        eta(1) = row * x - p.plant.x3d;
        for i = 2:r
            row    = row * sm.A;
            eta(i) = row * x;
        end
        % row is now Cy*A^(r-1)
        LgLfr1y = row * sm.B;
        Lfry    = row * sm.A * x;
        y       = eta(1);

    case 'pendulum'
        pv = p.plant.pv;
        [thd2, thd3, thd4f, thd4g] = ch5_pend_theta_lie(x, pv);

        y   = x(1:2) - p.plant.thetad(:);
        eta = [y; x(5:6); thd2; thd3];

        Lfry    = thd4f;
        LgLfr1y = thd4g;

    otherwise
        error('ch5_io_lin:unknownSystem', 'Unknown p.system "%s".', p.system);
end

rc = rcond(LgLfr1y);
if ~isfinite(rc) || rc < 1e-12
    % Neither plant should ever get here -- the springmass decoupling term is a
    % nonzero constant and the pendulum's is D^-1 k/Jm with D uniformly
    % positive definite -- so reaching this branch means the model changed, not
    % that the state is unusual. Fall back rather than emit NaN, and let
    % io.rcond say so.
    Ainv = pinv(LgLfr1y);
else
    Ainv = inv(LgLfr1y);
end

u_ff = -Ainv * Lfry;

io = struct('y', y, 'eta', eta, 'Lfry', Lfry, 'LgLfr1y', LgLfr1y, ...
            'Ainv', Ainv, 'u_ff', u_ff, 'rcond', rc);

end
