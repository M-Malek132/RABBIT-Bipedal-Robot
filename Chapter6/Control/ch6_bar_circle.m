function G = ch6_bar_circle(r, c, rho, sense, dc, ddc, drho, ddrho)
%CH6_BAR_CIRCLE  Signed distance to a (possibly moving) circle, to 2nd order.
%
%   G = ch6_bar_circle(r, c, rho, sense)
%   G = ch6_bar_circle(r, c, rho, sense, dc, ddc, drho, ddrho)
%
% Five of the six position constraints in Chapter 6 are "stay inside" or "stay
% outside" a circle -- O1 and O2 for footstep placement (6.9), the keep-out
% disc for an overhead obstacle (6.5), O3 and O4 for step width (6.20). They
% differ only in centre, radius and sense, so they are ONE function here and
% the individual barriers are three lines each.
%
%       outside :  g = |r - c| - rho >= 0
%       inside  :  g = rho - |r - c| >= 0
%
% Section 6.2 makes c and rho depend on time (the stones move), which is the
% reason this returns time partials as well as spatial ones. A static circle
% just leaves them zero.
%
% ------------------------------------------------------------------ the algebra
% With v = r - c, d = |v|, phi = v'/d, and Hd = (I - v v'/d^2)/d = d2d/dr2,
%
%       dg/dr    =  s phi
%       d2g/dr2  =  s Hd
%       dg/dt    =  s ( -phi cdot - rhodot )
%       d2g/dtdr =  s ( -cdot' Hd )
%       d2g/dt2  =  s ( cdot' Hd cdot - phi cddot - rhoddot )
%
% with s = +1 outside, -1 inside. The dg/dt row follows from dd/dt = -phi cdot
% at fixed r; the d2g/dt2 row from differentiating that once more, which is
% where the cdot' Hd cdot term comes from -- it is the curvature of the circle
% seen by a moving centre and it is the one term that is easy to drop. It is
% NOT small: for the sinusoidal stones of Fig. 6.10 it is the same order as the
% rest of d2g/dt2. ch6_test_barrier finite-differences all five against the
% definition rather than checking any of them by inspection.
%
% -------------------------------------------------------------- at the centre
% phi and Hd are undefined at r = c. That is not a numerical fragility to guard
% against with a fudge -- for an "inside" constraint the centre is the SAFEST
% point in the set, so a barrier row that has lost its gradient there is
% harmless, while for an "outside" constraint the centre is deep inside the
% forbidden region and the run is already over. Both are reported through
% G.singular so the caller decides; ch6_cbf_row drops the row and
% ch6_ctrl_cbf_clf_qp counts it.
%
% Inputs
%   r     : 2x1 point
%   c     : 2x1 centre
%   rho   : scalar radius
%   sense : 'inside' | 'outside'
%   dc, ddc     : 2x1 centre velocity / acceleration   (default 0)
%   drho, ddrho : scalar radius rate / acceleration    (default 0)
%
% Output
%   G : struct .g .dg (1x2) .Hg (2x2) .gt .gtr (1x2) .gtt .d .singular
%
% See also CH6_BAR_STONES, CH6_BAR_WIDTH, CH6_BAR_OBSTACLE, CH6_BAR_LIFT.

if nargin < 5 || isempty(dc),    dc    = zeros(2,1); end
if nargin < 6 || isempty(ddc),   ddc   = zeros(2,1); end
if nargin < 7 || isempty(drho),  drho  = 0;          end
if nargin < 8 || isempty(ddrho), ddrho = 0;          end

switch lower(sense)
    case 'outside', s =  1;
    case 'inside',  s = -1;
    otherwise
        error('ch6_bar_circle:sense', ...
              'sense must be inside|outside (got "%s").', sense);
end

v = r(:) - c(:);
d = norm(v);

G = struct('g', s*(d - rho), 'd', d, 'singular', false, ...
           'dg', zeros(1,2), 'Hg', zeros(2), ...
           'gt', 0, 'gtr', zeros(1,2), 'gtt', 0);

% Scale-relative, not absolute: a construction with rho ~ 1e-3 and one with
% rho ~ 1e3 do not deserve the same threshold.
if d <= 1e-12 * max(1, abs(rho))
    G.singular = true;
    return;
end

phi = v.' / d;                        % 1x2
Hd  = (eye(2) - (v*v.')/d^2) / d;     % 2x2

G.dg  = s * phi;
G.Hg  = s * Hd;
G.gt  = s * (-phi*dc(:) - drho);
G.gtr = s * (-(dc(:).' * Hd));
G.gtt = s * (dc(:).' * Hd * dc(:) - phi*ddc(:) - ddrho);

end
