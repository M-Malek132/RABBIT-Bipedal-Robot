function p = ch3_ub_defaults(p, mode, varargin)
%CH3_UB_DEFAULTS  Attach an upper-body wrench model to a Chapter-3 parameter struct.
%
%   p = ch3_ub_defaults(p, 'none')
%   p = ch3_ub_defaults(p, 'profile', 'fun', @(s) [Fx; Fz; M])
%   p = ch3_ub_defaults(p, 'bezier',  'beta', B, 'optimize', true, ...)
%
% WHAT IS MODELLED.  RABBIT's torso is one rigid 47 kg link (head, arms and
% trunk lumped, human-proportioned). The upper body is represented here as a
% WRENCH W = [Fx; Fz; M] acting on that torso:
%
%   Fx, Fz : world-frame force [N] (z up) applied at the "shoulder" point, a
%            distance l_sh from the hip along the torso axis;
%   M      : pure moment [N m], signed as the generalized force on q_t
%            (M > 0 pitches the torso toward larger q_t, i.e. forward).
%
% It enters the equations of motion as Q_ub = J_sh(q)' [Fx; Fz] + e_qt M.
% J_sh only has columns for px, pz, qt, so Q_ub(4:7) = 0 identically: the
% wrench acts DIRECTLY on x, z and q_t only and reaches the legs solely through
% the inertial coupling in M(q) and the stance-foot constraint.
%
% MODES
%   'none'    : no wrench (identical to not calling this function).
%   'profile' : W = fun(s), s in [0,1] the phase. A fixed, prescribed load --
%               a measured arm-swing reaction, a disturbance, a harness force.
%               Not a decision variable.
%   'bezier'  : W(s) = Bezier(beta, s), beta 3 x (deg+1). With optimize=true,
%               beta is appended to z and solved with the gait.
%
% PHYSICAL CAUTION.  Optimizing W as a free input is only meaningful when the
% source is genuinely EXTERNAL to the robot (crutch, handrail, harness, a
% therapist's hand). A robot's own arms/trunk are internal: they cannot change
% total momentum on their own, and a free W would act as a thruster. For
% internal upper-body motion model it as an extra link instead (see the
% accompanying notes). The bounds/sign constraints below are what keep an
% external W honest -- e.g. crutch Fz >= 0.
%
% Name/value options (defaults in brackets)
%   'l_sh'      shoulder distance from hip along torso [m]          [0.55]
%   'fun'       @(s) -> 3x1 wrench, for 'profile'                   [zeros]
%   'deg'       Bezier degree for 'bezier'                          [4]
%   'beta'      3 x (deg+1) initial/fixed coefficients              [zeros]
%   'optimize'  append beta to z                                    [false]
%   'lb','ub'   3x1 per-channel bounds on every coefficient [N;N;Nm]
%               (convex-hull property: bounds on the coefficients bound the
%               whole curve, so these are SUFFICIENT for W(s) in [lb, ub])
%                                                         [[-150;-150;-30], [150;150;30]]
%   'w'         weight on the wrench effort term added to the cost  [500]
%   'scale'     3x1 normalization of that term                      [[100;100;20]]
%   'periodic'  force W(0) = W(1) (continuous across steps)         [false]
%   'k_cone'    crutch-axis limit |Fx| <= k_cone * Fz, i.e. the force
%               points along a shaft at most atan(k_cone) from vertical
%               (0.4 -> 21.8 deg). [] = off. Only used with optimize.  [[]]
%
% See also CH3_UB_WRENCH, CH3_UB_COST, CH3_UB_INFO.

ip = inputParser;
ip.addParameter('l_sh', 0.55);
ip.addParameter('fun', @(s) zeros(3, 1));
ip.addParameter('deg', 4);
ip.addParameter('beta', []);
ip.addParameter('optimize', false);
ip.addParameter('lb', [-150; -150; -30]);
ip.addParameter('ub', [ 150;  150;  30]);
ip.addParameter('w', 500);
ip.addParameter('scale', [100; 100; 20]);
ip.addParameter('periodic', false);
ip.addParameter('k_cone', []);
ip.parse(varargin{:});
o = ip.Results;

ub = struct();
ub.mode     = lower(mode);
ub.l_sh     = o.l_sh;
ub.fun      = o.fun;
ub.optimize = logical(o.optimize);
ub.lb       = o.lb(:);
ub.ub       = o.ub(:);
ub.w        = o.w;
ub.scale    = o.scale(:);
ub.periodic = logical(o.periodic);
ub.k_cone   = o.k_cone;
if ~isempty(ub.k_cone) && ~(isscalar(ub.k_cone) && ub.k_cone >= 0 && isfinite(ub.k_cone))
    error('ch3_ub_defaults:k_cone', 'k_cone must be empty or a nonnegative scalar.');
end

if isempty(o.beta)
    ub.beta = zeros(3, o.deg + 1);
else
    if size(o.beta, 1) ~= 3
        error('ch3_ub_defaults:beta', 'beta must have 3 rows [Fx; Fz; M].');
    end
    ub.beta = o.beta;
end
if ~(o.l_sh >= 0 && o.l_sh <= 0.75)
    error('ch3_ub_defaults:l_sh', 'l_sh must lie on the 0.75 m torso (got %g).', o.l_sh);
end

p.ub = ub;
ch3_ub_info(p);            % validates mode / optimize combination
end
