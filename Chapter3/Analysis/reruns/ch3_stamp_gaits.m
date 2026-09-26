function T = ch3_stamp_gaits(files)
%CH3_STAMP_GAITS  Record the model in gait files that are orbits of it.
%
%   ch3_stamp_gaits()          every Results/ch3_gait_*.mat and the rerun gaits
%   ch3_stamp_gaits(files)     a cell array of paths
%   T = ch3_stamp_gaits(...)   one row per file
%
% Gaits solved before ch3_col_solve started storing out.model_sig carry no
% record of the dynamics they were solved on (ch3_model_check reports them
% 'unrecorded'). This is the one-off migration: each such file is re-verified
% on the dynamics now on the path, and ONLY a gait that is still an orbit --
% collocation equalities under 1e-4 (they sit at 1e-6 on a genuine gait and
% 0.3-0.7 on the five pre-2026-09-02 ones) AND a rollout within p.verify_tol of
% its nodes -- gets model_sig appended to its .mat. A gait that fails is left
% untouched and listed, so the stamp never certifies a stale file.
%
% Files already carrying a signature are only checked (match / mismatch).
%
% Output (also printed, and logged to Results/reruns/stamp_gaits.log)
%   T : struct array .file .status .max_ceq .verify_dev .stamped
%
% See also CH3_MODEL_SIGNATURE, CH3_MODEL_CHECK, CH3_COL_VERIFY.

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
if nargin < 1 || isempty(files)
    files = [list(fullfile(ROOT, 'Results', 'ch3_gait_*.mat')), ...
             list(fullfile(ROOT, 'Results', 'reruns', '*', 'gait_*.mat'))];
end
logf = fullfile(ROOT, 'Results', 'reruns', 'stamp_gaits.log');
if ~exist(fileparts(logf), 'dir'), mkdir(fileparts(logf)); end

EQ_TOL = 1e-4;
T = struct('file', {}, 'status', {}, 'max_ceq', {}, 'verify_dev', {}, 'stamped', {});
logln(logf, '=== ch3_stamp_gaits | %d files | %s', numel(files), datestr(now)); %#ok<TNOW1,DATST>

for k = 1:numel(files)
    f = files{k};
    row = struct('file', f, 'status', '', 'max_ceq', NaN, 'verify_dev', NaN, 'stamped', false);
    try
        S = load(f);
        [z, zname] = pick_z(S);
        if isempty(z) || ~isfield(S, 'p')
            row.status = 'no z / p';
            T(end+1) = row; logrow(logf, row, zname); continue; %#ok<AGROW>
        end
        R = ch3_model_check(S, 'quiet');
        if ~strcmp(R.status, 'unrecorded')
            row.status = R.status;               % already stamped: check only
            T(end+1) = row; logrow(logf, row, zname); continue; %#ok<AGROW>
        end
        p = ch3_upgrade_params(S.p);
        [~, ceq] = ch3_col_constraints(z, p);
        row.max_ceq = max(abs(ceq));
        V = ch3_col_verify(z, p, false);
        row.verify_dev = V.max_dev;
        if row.max_ceq <= EQ_TOL && V.ok
            model_sig = ch3_model_signature(p); %#ok<NASGU>
            save(f, 'model_sig', '-append');
            row.status  = 'orbit -> stamped';
            row.stamped = true;
        else
            row.status = 'NOT an orbit of this model -> left unstamped';
        end
    catch err
        row.status = ['error: ' err.message];
        zname = '?';
    end
    T(end+1) = row; %#ok<AGROW>
    logrow(logf, row, zname);
end
logln(logf, '=== STAMP_GAITS_DONE %d stamped of %d', sum([T.stamped]), numel(T));
fprintf('STAMP_GAITS_DONE\n');
end

% ---------------------------------------------------------------------------
function L = list(pattern)
F = dir(pattern);
L = arrayfun(@(d) fullfile(d.folder, d.name), F(:).', 'UniformOutput', false);
end

function [z, name] = pick_z(S)
z = []; name = '';
for c = {'z', 'z_opt', 'z_try'}
    if isfield(S, c{1}) && isnumeric(S.(c{1})) && ~isempty(S.(c{1}))
        z = S.(c{1})(:); name = c{1}; return;
    end
end
end

function logrow(logf, r, zname)
logln(logf, '%-70s %-5s max|ceq| %9.2e  verify %9.2e  %s', r.file, zname, ...
      r.max_ceq, r.verify_dev, r.status);
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
