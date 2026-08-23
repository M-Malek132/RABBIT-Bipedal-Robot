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
% ONE NESTED STRUCT DOES NOT TAKE THE DEFAULTS: p.limits.enable.  A gate the
% saved struct never mentioned is filled in as FALSE, not as whatever
% ch3_params says today, because its absence is evidence that the constraint
% was not enforced when that gait was solved rather than an opinion the result
% failed to record.  See merge_gates below for the eight-file incident that
% established this.
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

% Capture the gates AS SAVED before any merging: they are the only record of
% which constraints the solve actually saw, and merge() is about to bury them.
if isfield(p, 'limits') && isfield(p.limits, 'enable')
    saved_enable = p.limits.enable;
else
    saved_enable = struct();
end

p = merge(d, p);

if isfield(d, 'limits') && isfield(p, 'limits')
    p.limits = merge(d.limits, p.limits);
    if isfield(d.limits, 'enable') && isfield(p.limits, 'enable')
        p.limits.enable = merge_gates(d.limits.enable, saved_enable);
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
function out = merge_gates(defaults, saved)
% A GATE THE SAVED STRUCT NEVER MENTIONED WAS NOT ENFORCED.  So p.limits.enable
% is the one nested struct that does NOT take its missing fields from
% ch3_params: absent means OFF, and only an explicit saved true turns a
% constraint on.
%
% Everywhere else in this file the defaults are a reasonable stand-in for a
% field that did not exist yet -- a tolerance or a bound the result simply had
% no opinion about.  A gate is different.  Its absence is POSITIVE EVIDENCE
% about the solve: the constraint was not in the transcription when that gait
% was written, so nothing shaped the orbit to satisfy it, and the saved z is a
% minimizer of a problem that never contained the row.
%
% WHAT GOES WRONG OTHERWISE, MEASURED.  Results/ch3_gait_forward_lean_tall.mat
% was solved on 2026-07-28, when p.limits.enable had six fields and NEC3 did
% not exist.  b64160e added the six NIC/NEC gates defaulting to FALSE, so the
% merge still read that file correctly.  Then e7e101b flipped every default to
% TRUE, and from that commit on this function silently switched six
% constraints on underneath eight stored gaits.  The file then loaded claiming
% to enforce NEC3 while missing it by 0.92 (|Ix|/Iz = 0.463 against mu_s =
% 0.4) -- and it is the documented warm-start seed for ch3_lean_tall_march, so
% every warm run began infeasible in a direction the gait was never shaped
% for.  Nothing was written and nothing re-solved; a default in another file
% moved underneath them.
%
% ch3_col_verify does not catch this: it asks whether the nodes lie on a real
% trajectory, which they do (1.30e-05).  Being a real trajectory and being
% feasible for the current constraint set are different questions --
% ch3_col_check_limits is the one that asks the second.
%
% ENABLING A GATE IS ALWAYS AN EXPLICIT ACT, never something inherited from a
% default that moved.  ch3_impact_march says as much in its own header: it
% turns enable.impact on for its own copy of p precisely so the caller's stays
% off.  This makes the loader agree with that convention.
out = defaults;
fn = fieldnames(out);
for i = 1:numel(fn)
    out.(fn{i}) = false;
end
if isempty(fieldnames(saved)), return; end
sf = fieldnames(saved);
for i = 1:numel(sf)
    out.(sf{i}) = saved.(sf{i});
end
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
