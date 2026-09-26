function [f, g, aux] = ch4_control_affine(x, p, unc)
%CH4_CONTROL_AFFINE  Single-support dynamics of the TRUE (perturbed) model.
%
%   [f, g, aux] = ch4_control_affine(x, p)        uses p.uncertainty
%   [f, g, aux] = ch4_control_affine(x, p, unc)   uses the given perturbation
%   [f, g, aux] = ch4_control_affine(x, p, [])    NOMINAL model (no perturbation)
%
% Chapter 3 had one model. Chapter 4 has two, and the entire chapter is about
% the gap between them:
%
%   NOMINAL  ftil, gtil   the model the CONTROLLER was designed against
%   TRUE     f, g         the model the SIMULATION integrates
%
% This function produces either one. Passing unc = [] (or a perturbation that
% is identically trivial) short-circuits to ch3_control_affine, so the nominal
% model is not a re-derivation of Chapter 3 that could drift out of sync with
% it -- it IS Chapter 3, called directly.
%
% WHAT THE TWO PERTURBATIONS DO.
%
%   mass_scale s.  Scaling every link mass and inertia by s scales the kinetic
%   energy by s and the potential energy by s, so ALGEBRAICALLY M -> sM,
%   V -> sV, G -> sG (B unchanged, since actuators are not links). Rather than
%   apply that scaling to the nominal M/V/G in place, each scale the chapter
%   uses is rederived from scratch: rabbit_generate_case_dynamics reruns the
%   same Lagrangian trace construction that produced the nominal M.m/V.m/G.m,
%   with Mass_Properties_scaled(s) from the start, and writes a dedicated
%   M_<tag>/V_<tag>/G_<tag> per case. ch4_case_dynamics dispatches to them.
%   This does not take "linear in mass, so scaling commutes" on faith --
%   ch4_test_model compares the independently-rederived case against s*M(q)
%   etc. numerically.
%
%   Note what this does NOT do: it does not scale the torque limits. That is
%   the point of Case IV (s = 3), where a 3x robot must be driven by the same
%   actuators.
%
%   load_mass mL.  A point mass rigidly attached at the torso base. Its world
%   position is (px, -y) -- recall y is DOWN-positive -- so its Jacobian
%   Jp = [1 0 0 0 0 0 0; 0 -1 0 0 0 0 0] is CONSTANT. A constant Jacobian means
%   Jpdot = 0, so the load contributes
%
%       M += mL * Jp'Jp = mL * diag([1 1 0 0 0 0 0])
%       V += 0                        (no Coriolis: Jpdot = 0)
%       G += [0; -mL*g0; 0; 0; 0; 0; 0]
%
%   The gravity sign is fixed by the repo's convention U = m g0 z with z = -y,
%   which ch4_test_model checks against G(q) directly rather than assuming.
%
%   Unlike mass_scale this perturbation is NOT uniform, so it also changes the
%   impact map -- see ch4_impact.
%
%   links (with J_ref, b_visc, tau_coulomb, tau_bias).  A STRUCTURED model
%   error from ch4_uncertainty_set: each body's mass, COM and inertia on its
%   own (M, V, G from ch4_link_dynamics, exact for any link parameters),
%   reflected rotor inertia added to the actuated diagonal of M, and joint
%   friction and bias, which the robot adds to every command:
%
%       M qddot + V + G = B (u + tau_x(qdot)),    tau_x = ch4_joint_extra
%
%   tau_x depends on the state only, so it joins the DRIFT and the model stays
%   control-affine. A load_mass may ride on top; a mass_scale may not.
%
% Inputs
%   x   : 14x1 state [q; dq]
%   p   : parameter struct (uses p.nq, p.nu, p.g0)
%   unc : struct .mass_scale .load_mass; omitted -> p.uncertainty; [] -> none
%
% Outputs
%   f, g, aux : exactly as ch3_control_affine, but for the perturbed model.
%               aux additionally carries .unc, the perturbation actually used,
%               so downstream code can tell which model produced it.
%
% See also CH3_CONTROL_AFFINE, CH4_IMPACT, CH4_UNCERTAINTY.

if nargin < 3, unc = p.uncertainty; end

[s, mL] = unpack_unc(unc);
S = ch4_structured_part(unc);
if ~isempty(S) && s ~= 1
    error('ch4_control_affine:structuredScale', ...
          ['A structured model (uncertainty.links) cannot be combined with a ' ...
           'mass_scale; scale the link masses instead.']);
end

% --- nominal model: defer to Chapter 3 rather than reproduce it -----------
if s == 1 && mL == 0 && isempty(S)
    [f, g, aux] = ch3_control_affine(x, p);
    if nargout > 2
        aux.unc = struct('mass_scale', 1, 'load_mass', 0, 'structured', false);
    end
    return;
end

nq = p.nq;
nu = p.nu;

q  = x(1:nq);
dq = x(nq+1:2*nq);

% --- perturbed dynamics terms --------------------------------------------
tau_x = zeros(nu, 1);
if ~isempty(S)
    [M_mat, V_vec, G_vec] = ch4_link_dynamics(q, dq, S.links, p.g0);
    M_mat(4:7, 4:7) = M_mat(4:7, 4:7) + diag(S.J_ref);
    tau_x = ch4_joint_extra(dq, S);
elseif s == 1
    [M_mat, V_vec, G_vec] = ch3_mvg(q, dq, p);
else
    [M_mat, V_vec, G_vec] = ch4_case_dynamics(s, q, dq);
end
B_mat = input_matrix();

if mL ~= 0
    M_mat(1,1) = M_mat(1,1) + mL;
    M_mat(2,2) = M_mat(2,2) + mL;
    G_vec(2)   = G_vec(2)   - mL * p.g0;
end

% --- stance contact (geometry only: unaffected by the perturbation) -------
J      = J_st(q);
Jdotdq = Jdotdq_st(q, dq);
nc     = size(J, 1);

% --- one KKT factorization, five right-hand sides -------------------------
A = [M_mat, -J.'; ...
     J,     zeros(nc)];

rhs = [ [-V_vec - G_vec + B_mat * tau_x], B_mat ; ...
        [-Jdotdq],                         zeros(nc, nu) ];

sol = A \ rhs;

ddq_drift = sol(1:nq,     1);
ddq_in    = sol(1:nq,     2:end);
lam_drift = sol(nq+1:end, 1);
lam_in    = sol(nq+1:end, 2:end);

f = [dq;            ddq_drift];
g = [zeros(nq, nu); ddq_in   ];

if nargout > 2
    aux = struct('ddq_drift', ddq_drift, 'ddq_in', ddq_in, ...
                 'lam_drift', lam_drift, 'lam_in', lam_in, ...
                 'M', M_mat, 'Vv', V_vec, 'Gv', G_vec, ...
                 'J', J, 'Jdotdq', Jdotdq, 'tau_x', tau_x, ...
                 'unc', struct('mass_scale', s, 'load_mass', mL, ...
                               'structured', ~isempty(S)));
end

end

% ---------------------------------------------------------------------------
function [s, mL] = unpack_unc(unc)
%UNPACK_UNC  Normalize the perturbation spec, defaulting to "no perturbation".
if isempty(unc)
    s = 1; mL = 0;
    return;
end
if isfield(unc, 'mass_scale') && ~isempty(unc.mass_scale)
    s = unc.mass_scale;
else
    s = 1;
end
if isfield(unc, 'load_mass') && ~isempty(unc.load_mass)
    mL = unc.load_mass;
else
    mL = 0;
end
if ~(isscalar(s) && isfinite(s) && s > 0)
    error('ch4_control_affine:mass_scale', ...
          'uncertainty.mass_scale must be a positive finite scalar (got %s).', ...
          mat2str(s));
end
if ~(isscalar(mL) && isfinite(mL) && mL >= 0)
    error('ch4_control_affine:load_mass', ...
          'uncertainty.load_mass must be a nonnegative finite scalar (got %s).', ...
          mat2str(mL));
end
end
