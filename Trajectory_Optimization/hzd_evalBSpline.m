function [hd, dhd_dtheta] = hzd_evalBSpline(theta, CP, opt)
%HZD_EVALBSPLINE  Evaluate B-spline virtual constraint and its derivative.
%
%  Thin wrapper around your existing BSpline / BSpline_derivative functions.
%
%  Inputs:
%    theta  : current phase variable value (scalar)
%    CP     : (n+1)×4  control-point matrix
%    opt    : struct with fields  n_bs, p_bs, thetaStart, thetaEnd
%
%  Outputs:
%    hd           : 4×1  desired joint angles  h_d(s)
%    dhd_dtheta   : 4×1  d(h_d)/d(theta)
%
%  The B-spline is parameterised by the normalised phase
%      s = (theta - thetaStart) / (thetaEnd - thetaStart)  in [0,1]
%  so the chain rule gives:
%      dhd/dtheta = dhd/ds * ds/dtheta = dhd/ds / (thetaEnd - thetaStart)

n = opt.n_bs;
p = opt.p_bs;

% Normalised phase
s = (theta - opt.thetaStart) / (opt.thetaEnd - opt.thetaStart);
s = max(0.0, min(1.0, s));   % clamp to [0,1]

% Basis vector and its derivative  (both 1×(n+1))
N  = BSpline(n, p, s);           % from your Trajectory_Optimization/BSpline.m
dN = BSpline_derivative(n, p, s);% from your Trajectory_Optimization/BSpline_derivative.m

% Desired outputs:  CP is (n+1)×4,  N is 1×(n+1)
hd  = CP' * N';                  % 4×1

ds_dtheta   = 1.0 / (opt.thetaEnd - opt.thetaStart);
dhd_dtheta  = CP' * dN' * ds_dtheta;   % 4×1
end
