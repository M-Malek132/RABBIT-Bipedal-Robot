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
% Inputs
%   fname    : .mat file with a collocation solution (z and p). Default
%              Results/ch3_gait_upright.mat.
%   varargin : name/value overrides passed to ch4_params, applied AFTER the
%              gait fields are merged.
%
% Outputs
%   x0    : 14x1 start-of-step state of the periodic orbit
%   alpha : ny x n_ctrl virtual constraint coefficients
%   p     : Chapter-4 parameter struct consistent with that gait
%   meta  : struct .file .T .L_step .v_avg .source_p
%
% See also CH4_PARAMS, CH3_COL_UNPACK, CH4_MAIN.

if nargin < 1 || isempty(fname)
    here  = fileparts(mfilename('fullpath'));
    fname = fullfile(here, '..', 'Results', 'ch3_gait_upright.mat');
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

gait_fields = {'nq', 'nu', 'nx', 'iact', 'ny', 'H', ...
               'c_theta', 'theta_minus', 'theta_plus', ...
               'basis', 'bez_deg', 'n_ctrl', 'bsp_deg', 'g0'};

for i = 1:numel(gait_fields)
    f = gait_fields{i};
    if isfield(p_src, f)
        p.(f) = p_src.(f);
    end
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
% instead of the model error, and theta_hat runs away. Measured here at
% mass_scale 1 with no real uncertainty to find: against the 195 Nm posture
% gait a 65 Nm box gives theta_hat 2.6e5 and commanded torque 1.6e5 Nm; at
% 195 Nm the same run holds theta_hat at 5.2 and commands 198.8 Nm -- the
% same magnitude unconstrained 'l1' produces. The threshold is the gait's own
% peak and raising the box past it changes nothing, since it stops binding.
%
% The 65 Nm of Section 4.2.4 fits none of the gaits in Results/ (they run
% 195-465 Nm), so it is raised here rather than left to fail at run time. A
% caller that wants a smaller box still gets one -- the overrides below, and
% any later assignment to p.l1.u_max, win. That is what ch4_test_l1's
% saturation case relies on.
E = ch3_col_eval(S.z, ch3_upgrade_params(p_src));
p.gait_u_peak = max(max(abs(E.u(:))), max(abs(E.um(:))));

if p.l1.u_max < p.gait_u_peak
    warning('ch4_load_gait:l1Box', ...
            ['l1 torque box %.0f Nm is below this gait''s own peak %.1f Nm; ' ...
             'raising it. A box under the feedforward starves the inner QP ' ...
             'and the adaptation diverges.'], p.l1.u_max, p.gait_u_peak);
    p.l1.u_max = p.gait_u_peak;
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
                   'source_p', p_src);
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
