function G = ch6_bar_affine(r, a, a0, da0, dda0)
%CH6_BAR_AFFINE  A half-plane constraint  g = a0(t) + a'r >= 0, to 2nd order.
%
%   G = ch6_bar_affine(r, a, a0)
%   G = ch6_bar_affine(r, a, a0, da0, dda0)
%
% The one constraint in Chapter 6 that is not a circle: the low ceiling of
% (6.2), g_C = h_r - h_H, which is this with a = [0; -1] and a0 = h_r.
%
% It is worth keeping separate rather than modelling it as a circle of huge
% radius. A circle approximation of a plane has a curvature term Hg ~ 1/R that
% is small but not zero, and it enters gddot multiplied by |rdot|^2 -- at the
% 3.9 m/s the reference gait's foot reaches, a 100 m "flat" circle contributes
% 0.15 m/s^2 of spurious barrier acceleration. Exactly affine means exactly
% zero.
%
% Inputs
%   r   : 2x1 point
%   a   : 2x1 constant coefficient vector
%   a0  : scalar offset (may depend on time)
%   da0, dda0 : its time derivatives (default 0)
%
% Output
%   G : struct .g .dg (1x2) .Hg (2x2, zero) .gt .gtr (1x2, zero) .gtt .singular
%
% See also CH6_BAR_CIRCLE, CH6_BAR_OBSTACLE, CH6_BAR_LIFT.

if nargin < 4 || isempty(da0),  da0  = 0; end
if nargin < 5 || isempty(dda0), dda0 = 0; end

G = struct('g',   a0 + a(:).' * r(:), ...
           'd',   NaN, ...
           'singular', false, ...
           'dg',  a(:).', ...
           'Hg',  zeros(2), ...
           'gt',  da0, ...
           'gtr', zeros(1,2), ...
           'gtt', dda0);

end
