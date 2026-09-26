function S = ch4_structured_part(unc)
%CH4_STRUCTURED_PART  The structured model error in a p.uncertainty, or [].
%
%   S = ch4_structured_part(unc)
%
% p.uncertainty can describe the TRUE robot three ways: a uniform mass_scale,
% a load_mass at the hip, and -- from ch4_uncertainty_set -- per-link
% parameters (.links) with actuator terms. This returns the last, normalized,
% or [] when there is none, so ch4_control_affine and ch4_impact decide the
% same way which model the true robot is.
%
% Output (or []):
%   S.links       1x5 link parameters (ch4_link_params form)
%   S.J_ref       4x1 reflected rotor inertia on q1..q4 [kg m^2]
%   S.b_visc      4x1 viscous joint friction [N m s/rad]
%   S.tau_coulomb 4x1 Coulomb joint friction [N m]
%   S.tau_bias    4x1 constant torque bias [N m]
%   S.coulomb_vel scalar smoothing rate of the Coulomb term [rad/s]
%
% See also CH4_UNCERTAINTY_SET, CH4_JOINT_EXTRA, CH4_CONTROL_AFFINE.

S = [];
if isempty(unc) || ~isstruct(unc) || ~isfield(unc, 'links') || isempty(unc.links)
    return;
end
S = struct('links', unc.links, 'J_ref', zeros(4, 1), 'b_visc', zeros(4, 1), ...
           'tau_coulomb', zeros(4, 1), 'tau_bias', zeros(4, 1), ...
           'coulomb_vel', 0.05);
for f = {'J_ref', 'b_visc', 'tau_coulomb', 'tau_bias'}
    if isfield(unc, f{1}) && ~isempty(unc.(f{1}))
        v = unc.(f{1});
        if numel(v) == 1, v = repmat(v, 4, 1); end
        S.(f{1}) = v(:);
    end
end
if isfield(unc, 'coulomb_vel') && ~isempty(unc.coulomb_vel)
    S.coulomb_vel = unc.coulomb_vel(1);
end
if ~(S.coulomb_vel > 0)
    error('ch4_structured_part:coulombVel', ...
          'uncertainty.coulomb_vel must be positive (got %g).', S.coulomb_vel);
end
end
