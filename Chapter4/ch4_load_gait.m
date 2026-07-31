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
