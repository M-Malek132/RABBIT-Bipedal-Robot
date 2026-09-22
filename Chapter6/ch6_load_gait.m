function [x0, alpha, p, meta] = ch6_load_gait(fname, p)
%CH6_LOAD_GAIT  The nominal gait Chapter 6 runs on, checked against today's model.
%
%   [x0, alpha, p, meta] = ch6_load_gait()
%   [x0, alpha, p, meta] = ch6_load_gait(fname)
%   [x0, alpha, p, meta] = ch6_load_gait(fname, p)
%
% Chapter 6 designs no gait of its own for Sections 6.1-6.3: it bends ONE
% nominal Chapter-3 gait with barrier rows. Every Chapter-6 entry point starts
% here, so the choice of gait and the check that it is still a gait live in one
% place.
%
% THE DEFAULT IS Results/ch3_gait_posture_195.mat, the forward-lean 195 Nm gait
% solved on the 74 kg model -- the gait Chapters 3 and 4 are measured on.
% Chapter 6 was first built on Results/ch3_reference_gait.mat, which was solved
% before the 2026-09-02 regeneration of M/V/G (30 kg -> 74 kg). A gait file does
% not record which dynamics it was solved on, so that file still loads, unpacks
% and simulates -- as a reference that is no longer an orbit of the robot.
% meta.orbit re-evaluates the collocation residuals on the dynamics currently on
% the path, and this function refuses a gait that fails them unless the caller
% passes p.allow_stale_gait = true (only the regression that demonstrates the
% failure does).
%
% THE GAIT'S OWN PARAMETRIZATION COMES ALONG. alpha means nothing without the
% basis, degree and phase limits that produced it; loading it under ch6_params'
% defaults would silently describe a different curve. The fields that enter
% ch3_outputs / ch3_phase / ch3_yd are copied from the file's p into the
% returned p; every control knob stays Chapter 6's.
%
% Inputs
%   fname : .mat with a collocation solution (z or z_opt, and p). '' or missing
%           -> Results/ch3_gait_posture_195.mat
%   p     : a ch6_params struct to merge into (default ch6_params())
%
% Outputs
%   x0    : 14x1 start-of-step state of the orbit
%   alpha : ny x n_ctrl virtual-constraint coefficients
%   p     : p with the gait's parametrization merged in
%   meta  : .file .T .L_step .v_avg .u_peak .source_p .z .orbit
%           .orbit = .defect .periodicity .eta_post .tol .ok
%
% See also CH4_LOAD_GAIT, CH3_COL_EVAL, CH6_PARAMS, CH6_MAIN.

root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(fname)
    fname = fullfile(root, 'Results', 'ch3_gait_posture_195.mat');
end
if nargin < 2 || isempty(p), p = ch6_params(); end

if ~exist(fname, 'file')
    error('ch6_load_gait:notFound', ...
          'Gait file "%s" not found. Solve one with ch3_main.', fname);
end
S = load(fname);
if isfield(S, 'z'),         z = S.z;
elseif isfield(S, 'z_opt'), z = S.z_opt;
else
    error('ch6_load_gait:contents', '"%s" has no z or z_opt.', fname);
end
if ~isfield(S, 'p')
    error('ch6_load_gait:contents', '"%s" has no parameter struct p.', fname);
end
p_src = ch3_upgrade_params(S.p);

gait_fields = {'nq', 'nu', 'nx', 'iact', 'ny', 'H', ...
               'c_theta', 'theta_minus', 'theta_plus', ...
               'basis', 'bez_deg', 'n_ctrl', 'bsp_deg', 'g0'};
for i = 1:numel(gait_fields)
    f = gait_fields{i};
    if isfield(p_src, f), p.(f) = p_src.(f); end
end

% --- is it still an orbit of the dynamics on the path? ---------------------
E = ch3_col_eval(z, p_src);
orbit = struct('defect',      max(abs(E.defect(:))), ...
               'periodicity', norm(E.x_next(2:end) - E.X(2:end, 1)), ...
               'eta_post',    norm(E.eta_post), ...
               'tol',         p.verify_tol);
orbit.ok = max([orbit.defect, orbit.periodicity, orbit.eta_post]) <= orbit.tol;

if ~orbit.ok && ~(isfield(p, 'allow_stale_gait') && p.allow_stale_gait)
    error('ch6_load_gait:notAnOrbit', ...
          ['"%s" is not a periodic orbit of the dynamics currently on the ' ...
           'path (defect %.2e, periodicity %.2e, ||eta+|| %.2e; tol %.0e). ' ...
           'It was most likely solved before M/V/G were regenerated.'], ...
          fname, orbit.defect, orbit.periodicity, orbit.eta_post, orbit.tol);
end

[X, T, alpha] = ch3_col_unpack(z, p_src);
x0 = X(:, 1);

if nargout > 3
    meta = struct('file', fname, 'T', T, 'L_step', E.L_step, ...
                  'v_avg', E.L_step / T, ...
                  'u_peak', max(max(abs(E.u(:))), max(abs(E.um(:)))), ...
                  'source_p', p_src, 'z', z, 'orbit', orbit);
end
end
