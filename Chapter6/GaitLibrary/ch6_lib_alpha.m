function [alpha, info] = ch6_lib_alpha(lib, L_des)
%CH6_LIB_ALPHA  The gait library map  L_step -> alpha.  (6.22)-(6.24)
%
%   [alpha, info] = ch6_lib_alpha(lib, L_des)
%
% A finite set of gaits becomes a CONTINUUM by interpolating their Bezier
% coefficients linearly in step length:
%
%       zeta(L) = (L - L_i) / (L_(i+1) - L_i),      L_i <= L <= L_(i+1)   (6.22)
%       alpha(L) = (1 - zeta) alpha_i + zeta alpha_(i+1)                  (6.23)
%       A = { alpha(L) : L_1 <= L <= L_m }                                (6.24)
%
% and past the ends, linear EXTRAPOLATION off the last pair -- which is what the
% thesis does above 0.72 m and what Remark 6.6 then blames for part of the
% residual failure rate in Table 6.1. info.extrapolated says when it happened
% and info.extrap_frac by how much (in units of the end interval), so that
% attribution can be checked in a run rather than assumed.
%
% ================================================ WHY INTERPOLATION IS LEGITIMATE
% Interpolating the COEFFICIENTS is not the same as interpolating the gaits, and
% it is worth being clear which is claimed. A Bezier curve is linear in its
% coefficients, so alpha(L) generates exactly the pointwise convex combination
% of the two neighbouring output profiles -- that part is exact.
%
% What is NOT exact is that the resulting profile is a periodic orbit, or that
% it steps L. Nothing propagates the dynamics here. The thesis is explicit:
% "during steady-state, there is no theoretical guarantee that the interpolated
% gait results in exactly the desired step length", and none of the Table 3.1
% constraints are guaranteed during the transient after a switch. That gap is
% the entire reason Section 6.4 pairs the library with a CBF, and Remark 6.2
% names the typical failure mode -- foot scuffing during gait switches -- which
% is precisely what circle O2 rules out.
%
% So: the library provides a good GUESS at the next gait, and the barrier
% provides the guarantee. Neither does the other's job.
%
% Inputs
%   lib   : struct from ch6_lib_load -- .L (1 x m, ascending) and
%           .alpha (ny x n_ctrl x m)
%   L_des : desired step length [m]
%
% Output
%   alpha : ny x n_ctrl
%   info  : struct .i (lower bracket index) .zeta .extrapolated .extrap_frac
%           .L_lo .L_hi
%
% See also CH6_LIB_BUILD, CH6_LIB_LOAD, CH6_SIMULATE.

L = lib.L(:).';
m = numel(L);

if m == 0
    error('ch6_lib_alpha:empty', 'The gait library has no gaits.');
end
if m == 1
    alpha = lib.alpha(:,:,1);
    info  = struct('i', 1, 'zeta', 0, 'extrapolated', false, ...
                   'extrap_frac', 0, 'L_lo', L(1), 'L_hi', L(1));
    return;
end
if any(diff(L) <= 0)
    error('ch6_lib_alpha:order', ...
          'lib.L must be strictly ascending (got %s).', mat2str(L, 4));
end

% Bracket. The clamp to [1, m-1] is what turns interpolation into linear
% extrapolation outside the grid: zeta simply leaves [0,1] and (6.23) is
% evaluated unchanged.
i = find(L <= L_des, 1, 'last');
if isempty(i), i = 1; end
i = min(max(i, 1), m-1);

zeta = (L_des - L(i)) / (L(i+1) - L(i));                                 % (6.22)

alpha = (1 - zeta) * lib.alpha(:,:,i) + zeta * lib.alpha(:,:,i+1);       % (6.23)

extrap = (L_des < L(1)) || (L_des > L(m));
if L_des < L(1)
    frac = (L(1) - L_des) / (L(2) - L(1));
elseif L_des > L(m)
    frac = (L_des - L(m)) / (L(m) - L(m-1));
else
    frac = 0;
end

info = struct('i', i, 'zeta', zeta, 'extrapolated', extrap, ...
              'extrap_frac', frac, 'L_lo', L(i), 'L_hi', L(i+1));

end
