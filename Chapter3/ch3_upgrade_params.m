function p = ch3_upgrade_params(p)
%CH3_UPGRADE_PARAMS  Fill in fields a saved parameter struct predates.
%
%   p = ch3_upgrade_params(p)
%
% Merges current ch3_params defaults underneath the saved struct: saved
% fields win, missing fields are filled in. p.limits.enable is the exception
% -- a gate the saved struct never mentioned is filled in as FALSE, not
% today's default, since its absence means the constraint wasn't enforced
% when that gait was solved (see merge_gates). p.checkpoint_file is dropped
% rather than carried forward, since it names a scratch path on the machine
% that ran the solve, not a property of the gait.
%
% Call it at the top of anything that consumes a loaded p.
%
% See also CH3_PARAMS, CH3_COL_SOLVE.

d = ch3_params();

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
% Clear a checkpoint path this machine cannot honour ('' = no checkpointing).
if ~ischar(f) && ~isstring(f), f = ''; return; end
f = char(f);
if isempty(f), return; end

bad = false;
if ~ispc && (contains(f, '\') || ~isempty(regexp(f, '^[A-Za-z]:', 'once')))
    bad = true;             % Windows path seen from POSIX
end
if ~bad
    d = fileparts(f);
    if ~isempty(d) && ~isfolder(d), bad = true; end   % stale scratch dir
end
if bad, f = ''; end
end

% ---------------------------------------------------------------------------
function out = merge_gates(defaults, saved)
% Absent gate = not enforced, so it stays FALSE; only an explicit saved true
% turns a constraint on. See file header.
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
