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
% Call it at the top of anything that consumes a loaded p.
%
% See also CH3_PARAMS.

d = ch3_params();

p = merge(d, p);

if isfield(d, 'limits') && isfield(p, 'limits')
    p.limits = merge(d.limits, p.limits);
    if isfield(d.limits, 'enable') && isfield(p.limits, 'enable')
        p.limits.enable = merge(d.limits.enable, p.limits.enable);
    end
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
