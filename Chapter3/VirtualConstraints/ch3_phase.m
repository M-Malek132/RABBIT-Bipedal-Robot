function [s, ds_dq, theta, sdot] = ch3_phase(x, p)
%CH3_PHASE  Stage 2: the phase variable that replaces time.
%
%   [s, ds_dq, theta, sdot] = ch3_phase(x, p)
%
% Chapter 3 indexes the desired output profiles by a PHASE VARIABLE rather
% than by time:
%
%       theta = q_T + q_LA^st                     (absolute stance leg angle)
%       s     = (theta - theta_init) / (theta_final - theta_init)
%
% so s runs 0 -> 1 across a step.  This is what makes the gait SELF-PACED: the
% robot's own leg angle tells it how far through the step it is, and the
% resulting controller is time-invariant.  Push the robot and it re-indexes
% instead of falling off a pre-planned timeline.
%
% For RABBIT the absolute stance leg angle is
%
%       theta = qt + q1 + q2/2 = c*q,   c = [0 0 1 1 0.5 0 0]
%
% LINEAR in q by construction (see the note in ch3_params).  Two consequences
% the rest of the pipeline leans on:
%   * ds/dq is a CONSTANT row -- no state dependence to differentiate, which
%     is why Lf^2 y in ch3_io_lin has only one curvature term;
%   * theta cannot wrap, so the optimizer never sees a phantom jump.
%
% CLAMPING, AND WHY IT HAS A MARGIN.  A step that fails to strike can drive s
% far outside [0,1], where the polynomial extrapolates without bound, so there
% is a clamp -- but it engages only outside [-margin, 1+margin], NOT at the
% endpoints.
%
% Clamping hard at exactly [0,1] is a trap.  s = 0 is the start of every step,
% and floating point puts it a hair on either side; on the low side a hard
% clamp zeroes ds/dq, which zeroes sdot, which silently drops the (dyd/ds)sdot
% term from ydot.  The symptom is a gait that appears to start off the zero
% dynamics surface even when it was constructed to start on it.  A Bezier
% polynomial evaluated a whisker outside [0,1] is perfectly well behaved, so
% the margin costs nothing and removes a discontinuity that would otherwise
% fire at the start of every single step.
%
% Inputs
%   x : 14x1 state (or 7x1 q -- only the configuration is used for s)
%   p : parameter struct (uses p.c_theta, p.theta_minus, p.theta_plus, p.nq)
%
% Outputs
%   s     : scalar phase
%   ds_dq : 1x7 gradient ds/dq  (zero if the clamp is active)
%   theta : scalar raw phase variable
%   sdot  : scalar ds/dt = ds_dq * dq   (NaN if x was only a configuration)
%
% See also CH3_YD, CH3_OUTPUTS.

nq = p.nq;
q  = x(1:nq);

c     = p.c_theta(:).';
dth   = p.theta_plus - p.theta_minus;
theta = c * q;
s     = (theta - p.theta_minus) / dth;

ds_dq = c / dth;

do_clamp = ~isfield(p,'clamp_phase') || p.clamp_phase;
if isfield(p, 'phase_margin'), margin = p.phase_margin; else, margin = 0.5; end
if do_clamp && (s < -margin || s > 1 + margin)
    s     = min(max(s, -margin), 1 + margin);
    ds_dq = zeros(1, nq);
end

if nargout > 3
    if numel(x) >= 2*nq
        sdot = ds_dq * x(nq+1:2*nq);
    else
        sdot = NaN;
    end
end

end
