function [x0, alpha, p, meta] = ch4_load_gait(fname, varargin)
%CH4_LOAD_GAIT  Load a Chapter-3 optimized gait and set up Chapter-4 params.
%
%   [x0, alpha, p] = ch4_load_gait()
%   [x0, alpha, p] = ch4_load_gait(fname)
%   [x0, alpha, p] = ch4_load_gait(fname, 'uncertainty.mass_scale', 1.5, ...)
%
% Chapter 4 does not design a gait. It takes the periodic walking gait that
% Chapter 3 produced FOR THE NOMINAL MODEL and asks what happens when the robot
% is not that model. So every Chapter-4 entry point starts here.
%
% WHY THE GAIT'S OWN PARAMETERS MUST COME ALONG. alpha is meaningless without
% the parametrization that produced it: the Bezier degree fixes its shape, and
% the phase limits theta_minus / theta_plus fix what s = 0 and s = 1 mean. Load
% alpha under a different p and the virtual constraints silently describe a
% different curve -- the outputs would still evaluate, they would just be the
% wrong outputs. This function therefore copies the gait-defining fields out of
% the saved parameter struct and leaves only the CONTROL knobs at their
% Chapter-4 defaults.
%
% The fields copied are exactly those that enter ch3_outputs / ch3_phase /
% ch3_yd. Everything else -- controller choice, uncertainty, robust bounds, L1
% gains, sample rate -- is Chapter 4's to set.
%
% A GAIT FILE DOES NOT RECORD WHICH DYNAMICS IT WAS SOLVED ON. M.m/V.m/G.m are
% global, so a gait saved before they were regenerated still loads, still
% unpacks and still simulates -- as a reference that is no longer an orbit of
% the robot. That happened: the 2026-09-02 regeneration (30 kg -> 74 kg,
% non-uniformly) left the old default ch3_gait_upright.mat at a collocation
% defect of order one, and Chapter 4 ran on it for ten days, its robust sweep
% collapsing while every table still printed "3 steps completed". So the gait
% is re-evaluated here against TODAY's dynamics and meta.orbit reports whether
% it is still a periodic orbit (see the check below); ch4_main refuses one
% that is not. Gaits solved since then DO record their model (model_sig,
% stored by ch3_col_solve; ch3_stamp_gaits adds it to older files that still
% verify), and meta.model is ch3_model_check's verdict on it.
%
% Inputs
%   fname    : .mat file with a collocation solution (z and p). Default
%              Results/ch3_gait_posture_195.mat, the forward-lean 195 Nm gait
%              solved on the current 74 kg model.
%   varargin : name/value overrides passed to ch4_params, applied AFTER the
%              gait fields are merged.
%
% Outputs
%   x0    : 14x1 start-of-step state of the periodic orbit
%   alpha : ny x n_ctrl virtual constraint coefficients
%   p     : Chapter-4 parameter struct consistent with that gait
%   meta  : struct .file .T .L_step .v_avg .source_p .orbit .model
%           .orbit = .defect .periodicity .eta_post .tol .ok, the gait's own
%           collocation residuals under the dynamics currently on disk
%           .model = ch3_model_check: 'match' | 'mismatch' | 'unrecorded'
%
% See also CH4_PARAMS, CH3_COL_UNPACK, CH3_COL_EVAL, CH4_MAIN.

if nargin < 1 || isempty(fname)
    here  = fileparts(mfilename('fullpath'));
    fname = fullfile(here, '..', 'Results', 'ch3_gait_posture_195.mat');
end

if ~exist(fname, 'file')
    error('ch4_load_gait:notFound', ...
          ['Gait file "%s" not found. Solve one with ch3_main, or pass a ' ...
           'path to an existing Results/ch3_*.mat.'], fname);
end

S = load(fname);

if ~isfield(S, 'z') || ~isfield(S, 'p')
    error('ch4_load_gait:contents', ...
          '"%s" must contain both z (decision vector) and p (its parameters).', ...
          fname);
end

p_src = S.p;

% --- Chapter-4 defaults, then the gait's own parametrization --------------
p = ch4_params();

% free_theta says whether z carries the phase endpoints (ch3_col_pack); a gait
% solved with them free is saved with its solved theta_minus / theta_plus in p.
gait_fields = {'nq', 'nu', 'nx', 'iact', 'ny', 'H', ...
               'c_theta', 'theta_minus', 'theta_plus', 'free_theta', ...
               'basis', 'bez_deg', 'n_ctrl', 'bsp_deg', 'g0'};

for i = 1:numel(gait_fields)
    f = gait_fields{i};
    if isfield(p_src, f)
        p.(f) = p_src.(f);
    end
end

% --- which dynamics was it solved on? -------------------------------------
% A file written since ch3_col_solve started storing model_sig says so, and a
% mismatch is conclusive. A file from before (posture_195 until
% ch3_stamp_gaits has run) is 'unrecorded', and the orbit check below is what
% decides for it.
model = ch3_model_check(S, 'quiet');
if strcmp(model.status, 'mismatch')
    warning('ch4_load_gait:modelMismatch', '"%s": %s', fname, model.msg);
end

% --- re-evaluate the gait on the dynamics on the path TODAY --------------
E = ch3_col_eval(S.z, ch3_upgrade_params(p_src));

% IS IT STILL AN ORBIT? ch3_col_eval recomputes the collocation defects, the
% periodicity through Delta (px excluded, as in ch3_col_constraints) and the
% post-impact outputs with whatever M/V/G are current. A gait solved on those
% dynamics has all three at solver tolerance -- posture_195 reads 7.6e-07,
% 7.6e-09 and 1.2e-06 -- while ch3_gait_upright, solved before the
% regeneration, reads 6.6e-02, 1.02 and 0.51. p.verify_tol, Chapter 3's own
% threshold for "not a real trajectory", separates the two by three orders of
% magnitude either way.
orbit = struct('defect',      max(abs(E.defect(:))), ...
               'periodicity', norm(E.x_next(2:end) - E.X(2:end, 1)), ...
               'eta_post',    norm(E.eta_post), ...
               'tol',         p.verify_tol);
orbit.ok = max([orbit.defect, orbit.periodicity, orbit.eta_post]) <= orbit.tol;

if ~orbit.ok
    warning('ch4_load_gait:notAnOrbit', ...
            ['"%s" is not a periodic orbit of the dynamics currently on the ' ...
             'path (defect %.2e, periodicity %.2e, ||eta+|| %.2e; tol %.0e). ' ...
             'It was most likely solved before M/V/G were regenerated; ' ...
             're-solve it in Chapter 3 before measuring controllers on it.'], ...
            fname, orbit.defect, orbit.periodicity, orbit.eta_post, orbit.tol);
end

% --- size the L1 torque box to the gait ----------------------------------
% p.gait_u_peak is the largest joint torque the gait itself uses. It is the
% floor under any saturation box: the CLF-QP inside 'l1_con' has to deliver
% the FEEDFORWARD before it can correct anything, so a box below this makes
% the inner QP infeasible at most samples rather than merely tight.
%
% WHY THAT MATTERS MORE FOR L1 THAN FOR THE ROBUST LAWS. A starved
% clfqp_con/rclfqp_con just tracks badly and stays bounded. L1 does not: the
% adaptation is driven by -G'Peps eta_tilde ||eta|| at Gamma = 1e4, so once
% the robot leaves the orbit the estimator chases the TRACKING FAILURE
% instead of the model error, and theta_hat runs away. Measured here (at the
% then-default eps = 0.35, as is every number in this note) at
% mass_scale 1 with no real uncertainty to find: against the 195 Nm posture
% gait a 65 Nm box gives theta_hat 2.6e5 and commanded torque 1.6e5 Nm; at
% 195 Nm the same run holds theta_hat at 5.2 and commands 198.8 Nm -- the
% same magnitude unconstrained 'l1' produces.
%
% THE FLOOR NEEDS HEADROOM, NOT JUST THE PEAK. At mass_scale 1 the box binds
% only on the feedforward. Under a real perturbation mu1 must carry feedback
% as well until the estimate catches up, and a box sitting exactly on the
% feedforward peak leaves it none. Measured at mass_scale 1.5, 3 steps: box at
% 1.00x the peak ran away to 7193 Nm and lost the robot in step 3; at 1.25x
% it walked all three at max||eta|| 7.7, Vend/Vmx 0.065; 1.5x and 2x changed
% nothing material (Vend/Vmx 0.036, 0.040). Mass scale 0.7 degrades to
% max||eta|| ~21 at EVERY box, unconstrained 'l1' included, so that case is
% not the box.
%
% The 65 Nm of Section 4.2.4 fits none of the gaits in Results/ (they run
% 195-465 Nm), so it is raised here rather than left to fail at run time. A
% caller that wants a smaller box still gets one -- the overrides below, and
% any later assignment to p.l1.u_max, win. That is what ch4_test_l1's
% saturation case relies on.
L1_HEADROOM = 1.25;

p.gait_u_peak = max(max(abs(E.u(:))), max(abs(E.um(:))));

if p.l1.u_max < L1_HEADROOM * p.gait_u_peak
    warning('ch4_load_gait:l1Box', ...
            ['l1 torque box %.0f Nm is below %.2fx this gait''s own peak ' ...
             '(%.1f Nm); raising it to %.0f Nm. A box without headroom over ' ...
             'the feedforward starves the inner QP and the adaptation ' ...
             'diverges.'], p.l1.u_max, L1_HEADROOM, p.gait_u_peak, ...
            L1_HEADROOM * p.gait_u_peak);
    p.l1.u_max = L1_HEADROOM * p.gait_u_peak;
end

% --- overrides last, so a caller can still change anything ---------------
if ~isempty(varargin)
    p = ch4_params_override(p, varargin);
end

[X, T, alpha] = ch3_col_unpack(S.z, p);
x0 = X(:, 1);

if nargout > 3
    foot0 = P_st(X(1:p.nq, 1));
    footN = P_sw(X(1:p.nq, end));
    L     = footN(1) - foot0(1);
    meta  = struct('file', fname, 'T', T, 'L_step', L, 'v_avg', L / T, ...
                   'source_p', p_src, 'orbit', orbit, 'model', model);
end

end

% ---------------------------------------------------------------------------
function p = ch4_params_override(p, nv)
%CH4_PARAMS_OVERRIDE  Apply name/value pairs, including dotted nested names.
for k = 1:2:numel(nv)
    p = set_field(p, nv{k}, nv{k+1});
end
end

function s = set_field(s, name, value)
parts = strsplit(name, '.');
if ~isfield(s, parts{1})
    error('ch4_load_gait:unknownField', 'Unknown parameter "%s".', name);
end
if numel(parts) == 1
    s.(parts{1}) = value;
else
    s.(parts{1}) = set_field(s.(parts{1}), strjoin(parts(2:end), '.'), value);
end
end
