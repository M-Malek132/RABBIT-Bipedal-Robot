function R = ch3_model_check(src, mode)
%CH3_MODEL_CHECK  Was this gait solved on the dynamics that are on the path now?
%
%   R = ch3_model_check(S)          S = load(gait_file)
%   R = ch3_model_check(p)          a parameter struct carrying p.model_sig
%   R = ch3_model_check(src, mode)  mode: 'quiet' | 'report' (default) |
%                                         'warn' | 'error'
%
% The signature a solve stored (ch3_col_solve: out.model_sig and
% out.p.model_sig; ch3_stamp_gaits for files that predate it) is looked for
% in S.model_sig, S.out.model_sig, S.solve_out.model_sig, S.p.model_sig, then
% src.model_sig, and compared with ch3_model_signature of the SAME model
% variant (its p.model_blend) evaluated now. Values, not digests, are compared
% (relative tolerance 1e-9), so the check survives last-bit differences
% between machines.
%
%   'match'       the file's model is the one on the path
%   'mismatch'    it is NOT: the gait was solved on other dynamics, and every
%                 number computed from it here is about a different robot
%   'unrecorded'  the file predates signatures; only re-verifying it on
%                 today's model (ch3_col_verify, or ch4_load_gait's orbit
%                 check) can say whether it is still an orbit
%
% mode 'warn' warns on mismatch and unrecorded; 'error' raises on mismatch
% (id ch3_model_check:mismatch) and warns on unrecorded; 'report' prints one
% line; 'quiet' prints nothing.
%
% Output
%   R : struct .status .rel_err .saved .current .msg
%
% See also CH3_MODEL_SIGNATURE, CH3_STAMP_GAITS, CH4_LOAD_GAIT.

if nargin < 2 || isempty(mode), mode = 'report'; end
TOL = 1e-9;

[sig, p_src] = find_signature(src);

R = struct('status', 'unrecorded', 'rel_err', NaN, 'saved', sig, ...
           'current', [], 'msg', '');

if isempty(sig)
    R.msg = ['the gait file does not record the dynamics it was solved on ' ...
             '(it predates ch3_model_signature). Re-verify it on today''s ' ...
             'model before trusting it; ch3_stamp_gaits records the ' ...
             'signature once it passes.'];
else
    pc = struct('model_blend', [], 'g0', 9.8062);
    if isfield(sig, 'blend'), pc.model_blend = sig.blend; end
    if ~isempty(p_src) && isfield(p_src, 'g0'), pc.g0 = p_src.g0; end
    cur = ch3_model_signature(pc);
    R.current = cur;
    if ~isfield(sig, 'version') || sig.version ~= cur.version || ...
            numel(sig.vals) ~= numel(cur.vals)
        R.msg = ['the stored signature has another layout (version) and ' ...
                 'cannot be compared; treat the gait as unrecorded.'];
    else
        R.rel_err = max(abs(sig.vals(:) - cur.vals(:))) / max(1, max(abs(cur.vals(:))));
        if R.rel_err <= TOL
            R.status = 'match';
            R.msg = sprintf('solved on the model on the path (%s, %.1f kg)', ...
                            cur.digest, cur.mass);
        else
            R.status = 'mismatch';
            R.msg = sprintf(['solved on OTHER dynamics: model %s (%.1f kg) in the ' ...
                             'file, %s (%.1f kg) on the path, relative ' ...
                             'difference %.2e. This gait is not an orbit of ' ...
                             'the robot being simulated.'], ...
                            sig.digest, sig.mass, cur.digest, cur.mass, R.rel_err);
        end
    end
end

switch lower(mode)
    case 'quiet'
    case 'report'
        fprintf(' MODEL   %-10s %s\n', R.status, R.msg);
    case 'warn'
        if ~strcmp(R.status, 'match')
            warning(['ch3_model_check:' R.status], '%s', R.msg);
        end
    case 'error'
        if strcmp(R.status, 'mismatch')
            error('ch3_model_check:mismatch', '%s', R.msg);
        elseif ~strcmp(R.status, 'match')
            warning('ch3_model_check:unrecorded', '%s', R.msg);
        end
    otherwise
        error('ch3_model_check:mode', ...
              'mode must be quiet | report | warn | error (got "%s").', mode);
end

end

% ---------------------------------------------------------------------------
function [sig, p_src] = find_signature(src)
sig = []; p_src = [];
if ~isstruct(src), return; end
if isfield(src, 'p') && isstruct(src.p), p_src = src.p; else, p_src = src; end
cands = {};
if isfield(src, 'model_sig'), cands{end+1} = src.model_sig; end
for f = {'out', 'solve_out'}
    if isfield(src, f{1}) && isstruct(src.(f{1})) && isfield(src.(f{1}), 'model_sig')
        cands{end+1} = src.(f{1}).model_sig; %#ok<AGROW>
    end
end
if isstruct(p_src) && isfield(p_src, 'model_sig')
    cands{end+1} = p_src.model_sig;
end
for k = 1:numel(cands)
    c = cands{k};
    if isstruct(c) && isfield(c, 'vals') && ~isempty(c.vals)
        sig = c;
        return;
    end
end
end
