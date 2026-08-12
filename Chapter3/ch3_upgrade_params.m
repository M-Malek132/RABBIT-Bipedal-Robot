function p = ch3_upgrade_params(p)
%CH3_UPGRADE_PARAMS  Fill in fields a saved parameter struct predates.
%
%   p = ch3_upgrade_params(p)
%
% Results are saved with the parameter struct that produced them, so a .mat
% written last week carries last week's field list.  Loading it and calling
% anything that reads a field added since fails with "Unrecognized field
% name", which looks like a bug in the analysis code rather than what it is --
% an old file.
%
% This merges the CURRENT defaults underneath the saved values: every field the
% saved struct has is kept exactly as saved (so a result is still analysed with
% the settings that produced it), and only genuinely missing fields are filled
% in from ch3_params.  Nested limit structs are merged the same way, since new
% limits get added there too.
%
% ONE FIELD IS DELIBERATELY NOT PRESERVED: p.checkpoint_file.  Everything else
% here describes the GAIT, and re-analysing a result under the settings that
% produced it is the whole point; checkpoint_file describes nothing but a
% scratch file on the machine that ran the solve, so carrying it forward has no
% upside and one sharp edge.  A .mat saved on Windows stores an absolute
% 'C:\Users\...\Results\ch3_*_ckpt.mat', and on macOS/Linux a backslash is a
% LEGAL FILENAME CHARACTER -- so ch3_col_solve's checkpoint save does not fail,
% it succeeds, creating one file whose entire Windows path is its name in
% whatever the working directory happened to be.  Two of those were committed to
% this repo before anyone noticed.  ch3_col_solve guards the save with a
% try/catch, but that only catches an UNWRITABLE path; this one is writable and
% simply wrong, so it has to be caught here instead.
%
% Call it at the top of anything that consumes a loaded p.
%
% See also CH3_PARAMS, CH3_COL_SOLVE.

d = ch3_params();

p = merge(d, p);

if isfield(d, 'limits') && isfield(p, 'limits')
    p.limits = merge(d.limits, p.limits);
    if isfield(d.limits, 'enable') && isfield(p.limits, 'enable')
        p.limits.enable = merge(d.limits.enable, p.limits.enable);
    end
end

p.checkpoint_file = usable_ckpt(p.checkpoint_file);

end

% ---------------------------------------------------------------------------
function f = usable_ckpt(f)
% Empty out a checkpoint path this machine cannot honour.  '' is the documented
% "no checkpointing" value, and ch3_main fills a fresh one in when it is empty,
% so clearing is always safe -- worst case the run simply goes uncheckpointed.
if ~ischar(f) && ~isstring(f), f = ''; return; end
f = char(f);
if isempty(f), return; end

bad = false;

% A Windows path seen from POSIX: drive letter, or backslash separators that
% this filesystem would swallow into the filename rather than reject.
if ~ispc && (contains(f, '\') || ~isempty(regexp(f, '^[A-Za-z]:', 'once')))
    bad = true;
end

% Otherwise the usual stale-scratch case: the directory is simply gone.
if ~bad
    d = fileparts(f);
    if ~isempty(d) && ~isfolder(d), bad = true; end
end

if bad, f = ''; end
end

% ---------------------------------------------------------------------------
function out = merge(defaults, given)
% Start from the defaults, then let every field the caller actually has win.
out = defaults;
fn = fieldnames(given);
for i = 1:numel(fn)
    out.(fn{i}) = given.(fn{i});
end
end
