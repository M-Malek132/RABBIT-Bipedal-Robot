function B = ch4_delta_bounds_prior(X, alpha, p, U, opts)
%CH4_DELTA_BOUNDS_PRIOR  D1, D2 as maxima over a PRIOR set of model errors.
%
%   B = ch4_delta_bounds_prior(X, alpha, p, U)
%   B = ch4_delta_bounds_prior(X, alpha, p, U, opts)
%
% ch4_delta_bounds measures Delta1, Delta2 under the very perturbation the
% robust law is then run against, so the law is told the size of the
% disturbance it will face -- an oracle, not a design. Here the bounds come
% from a prior set instead (ch4_uncertainty_set): the maximum of ||Delta1||
% and ||Delta2|| over draws from the set, over the states X (the gait's own
% nodes) and, optionally, a jittered cloud around them, times a safety factor.
% The robust law run with these numbers knows the SET, never the draw it
% meets -- which is what a deployed controller could know.
%
% The maximum over a finite sample of the set is an estimate of the supremum,
% not a guarantee; the draws' own maxima are returned so the spread (and how
% close the sample came to its largest value) can be read off. The two
% extreme corners -- every mass, inertia and COM offset at the top of its
% range, then at the bottom, actuator terms at their largest -- are always
% included.
%
% Inputs
%   X     : 14 x nt states (e.g. the gait's collocation nodes)
%   alpha : ny x n_ctrl
%   p     : parameter struct
%   U     : 1 x n draws from ch4_uncertainty_set
%   opts  : .safety (1.2, as the chapter's own bounds) .jitter (0.05)
%           .n_jitter (0) .spec (the set U was drawn from, for the corners)
%
% Output
%   B : struct .delta1_max .delta2_max (with safety) .n1_max .n2_max
%       .n1 .n2 (per draw) .n2_scalar_max .feasible (delta2_max < 1)
%
% See also CH4_UNCERTAINTY_SET, CH4_UNCERTAINTY, CH4_DELTA_BOUNDS.

if nargin < 5, opts = struct(); end
d = struct('safety', 1.2, 'jitter', 0.05, 'n_jitter', 0, 'spec', []);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = d.(f{i}); end
end
spec = opts.spec;
if isempty(spec), spec = ch4_uncertainty_set(); end

U = [U, corners(spec)];
rs = RandStream('mt19937ar', 'Seed', 7);          % reproducible jitter cloud
nt = size(X, 2);
n1 = zeros(1, numel(U)); n2 = n1; n2s = n1;
for k = 1:numel(U)
    for j = 1:nt
        for r = 0:opts.n_jitter
            x = X(:, j);
            if r > 0, x = x + opts.jitter * randn(rs, size(x)); end
            D = ch4_uncertainty(x, alpha, p, U(k));
            if ~isfinite(D.n1) || ~isfinite(D.n2), continue; end
            n1(k)  = max(n1(k),  D.n1);
            n2(k)  = max(n2(k),  D.n2);
            n2s(k) = max(n2s(k), D.n2_scalar);
        end
    end
end

B = struct('n1', n1, 'n2', n2, 'n1_max', max(n1), 'n2_max', max(n2), ...
           'n2_scalar_max', max(n2s), 'safety', opts.safety, ...
           'delta1_max', opts.safety * max(n1), ...
           'delta2_max', opts.safety * max(n2), 'n_draws', numel(U), ...
           'n_states', nt * (1 + opts.n_jitter));
B.feasible = B.delta2_max < 1;
end

% ---------------------------------------------------------------------------
function C = corners(spec)
%CORNERS  Every link parameter at the top of its range, then at the bottom.
L0 = ch4_link_params();
C  = ch4_uncertainty_set(2, 1, spec);                 % for the field layout
for c = 1:2
    sgn = 3 - 2*c;                                     % +1, then -1
    L = L0;
    for i = 1:5
        L(i).m = L0(i).m * (1 + sgn * spec.mass_rel);
        L(i).b = L0(i).b + sign(L0(i).b) * sgn * spec.com_abs;
        L(i).J = L0(i).J * (1 + sgn * spec.inertia_rel);
    end
    C(c).links       = L;
    C(c).J_ref       = spec.J_ref       * ones(4, 1);
    C(c).b_visc      = spec.b_visc      * ones(4, 1);
    C(c).tau_coulomb = spec.tau_coulomb * ones(4, 1);
    C(c).tau_bias    = sgn * spec.tau_bias * ones(4, 1);
end
end
