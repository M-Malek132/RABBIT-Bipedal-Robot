function [U, spec] = ch4_uncertainty_set(n, seed, spec)
%CH4_UNCERTAINTY_SET  A STRUCTURED prior set of model errors, and draws from it.
%
%   spec     = ch4_uncertainty_set()                the default set, no draws
%   [U,spec] = ch4_uncertainty_set(n, seed)         n draws, reproducibly
%   [U,spec] = ch4_uncertainty_set(n, seed, spec)   from a modified set
%
% WHY. Chapter 4's robustness test is a uniform mass scale -- which, as the
% chapter itself derives, is exactly an input-gain error plus a proportional
% drift error, with the zero dynamics unchanged: about the easiest uncertainty
% an adaptive law can face -- and a hip load, where the claims shrink. And its
% robust bounds D1, D2 are measured FROM the disturbance being tested, so the
% robust law is told the answer. A structured prior set fixes both: model
% errors a real robot has, set BEFORE any run, with D1 and D2 computed as
% maxima over the set (ch4_delta_bounds_prior) and the controllers then tested
% on fresh draws from it (ch4_structured_study).
%
% THE SET (every range is an ASSUMPTION of this study, not an identified
% value; change spec to change it):
%
%   links        each of the five bodies independently (spec.symmetric =
%                false; true perturbs both legs alike):
%     mass_rel     mass  x (1 + U[-r, r])                    r = 0.10
%     com_abs      COM shifted along the link by U[-d, d] m  d = 0.02
%     inertia_rel  planar inertia x (1 + U[-r, r])           r = 0.20
%   actuators, per joint (q1..q4):
%     J_ref        reflected rotor inertia U[0, J] kg m^2    J = 0.20
%                  RABBIT drives each joint through a 50:1 harmonic drive
%                  (Chevallereau et al. 2003, Table I); a rotor of inertia
%                  J_r appears at the joint as 2500 J_r, and the model has none
%     b_visc       viscous friction U[0, b] N m s/rad        b = 1.0
%     tau_coulomb  Coulomb friction U[0, c] N m, smoothed    c = 2.0
%                  as c tanh(qdot / coulomb_vel), coulomb_vel = 0.05 rad/s
%     tau_bias     a torque bias U[-t, t] N m                t = 2.0
%   implementation (the simulator's, not the model's):
%     noise_q      measurement noise on q, std [rad]         1e-3
%     noise_dq     ... on qdot, std [rad/s]                  1e-2
%     delay        actuation delay in control samples        1
%
% Output
%   spec : the set actually used
%   U    : 1 x n struct array of p.uncertainty entries, each with
%          .mass_scale = 1 .load_mass = 0 .links (1x5, ch4_link_params form)
%          .J_ref .b_visc .tau_coulomb .tau_bias (4x1) .coulomb_vel
%          .noise (struct .q .dq .seed) .delay .asymmetric
%
% See also CH4_LINK_PARAMS, CH4_LINK_DYNAMICS, CH4_CONTROL_AFFINE,
%          CH4_DELTA_BOUNDS_PRIOR, CH4_STRUCTURED_STUDY.

d = struct('mass_rel', 0.10, 'com_abs', 0.02, 'inertia_rel', 0.20, ...
           'J_ref', 0.20, 'b_visc', 1.0, 'tau_coulomb', 2.0, ...
           'coulomb_vel', 0.05, 'tau_bias', 2.0, ...
           'noise_q', 1e-3, 'noise_dq', 1e-2, 'delay', 1, ...
           'symmetric', false);
if nargin < 3 || isempty(spec)
    spec = d;
else
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(spec, f{i}) || isempty(spec.(f{i})), spec.(f{i}) = d.(f{i}); end
    end
end

if nargin < 1 || isempty(n)
    U = spec;                 % one-output form: the set itself
    return;
end
if nargin < 2 || isempty(seed), seed = 1; end

rs = RandStream('mt19937ar', 'Seed', seed);
L0 = ch4_link_params();
U  = repmat(empty_entry(), 1, n);
for k = 1:n
    L = L0;
    draws = rand(rs, 5, 3) * 2 - 1;                    % U[-1, 1] per link
    if spec.symmetric
        draws(4, :) = draws(2, :);                     % swing thigh = stance
        draws(5, :) = draws(3, :);                     % swing shank = stance
    end
    for i = 1:5
        L(i).m = L0(i).m * (1 + spec.mass_rel    * draws(i, 1));
        L(i).b = L0(i).b + sign(L0(i).b) * spec.com_abs * draws(i, 2);
        L(i).J = L0(i).J * (1 + spec.inertia_rel * draws(i, 3));
    end
    jd = rand(rs, 4, 4);
    if spec.symmetric
        jd(3:4, :) = jd(1:2, :);                       % swing joints = stance
    end
    e = empty_entry();
    e.links       = L;
    e.J_ref       = spec.J_ref       * jd(:, 1);
    e.b_visc      = spec.b_visc      * jd(:, 2);
    e.tau_coulomb = spec.tau_coulomb * jd(:, 3);
    e.tau_bias    = spec.tau_bias    * (2 * jd(:, 4) - 1);
    e.coulomb_vel = spec.coulomb_vel;
    e.noise       = struct('q', spec.noise_q, 'dq', spec.noise_dq, ...
                           'seed', seed * 1000 + k);
    e.delay       = spec.delay;
    e.asymmetric  = ~spec.symmetric;
    U(k) = e;
end
end

% ---------------------------------------------------------------------------
function e = empty_entry()
e = struct('mass_scale', 1, 'load_mass', 0, 'links', ch4_link_params(), ...
           'J_ref', zeros(4, 1), 'b_visc', zeros(4, 1), ...
           'tau_coulomb', zeros(4, 1), 'tau_bias', zeros(4, 1), ...
           'coulomb_vel', 0.05, ...
           'noise', struct('q', 0, 'dq', 0, 'seed', 0), ...
           'delay', 0, 'asymmetric', false);
end
