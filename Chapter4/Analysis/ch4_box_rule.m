function rule = ch4_box_rule(p)
%CH4_BOX_RULE  Which torque-box rule the sweeps use: 'rating' or 'thesis'.
%
%   rule = ch4_box_rule(p)
%
% p.box.rule, lower-cased and checked (see ch4_params). A parameter struct
% saved before p.box existed -- an old ch4_result .mat -- has no field and
% resolves to 'thesis', the rule that run actually used, so re-analysing it
% reproduces its boxes rather than silently swapping in the rating.
%
% See also CH4_PARAMS, CH4_COMPARE_CONTROLLERS, CH4_LOAD_STUDY.

if ~isfield(p, 'box') || ~isfield(p.box, 'rule') || isempty(p.box.rule)
    rule = 'thesis';
    return;
end

rule = lower(p.box.rule);
if ~any(strcmp(rule, {'rating', 'thesis'}))
    error('ch4_box_rule:rule', ...
          'Unknown p.box.rule "%s" (expected rating|thesis).', p.box.rule);
end
if strcmp(rule, 'rating') && ~(isfield(p.box, 'rating') && isscalar(p.box.rating) ...
                               && isfinite(p.box.rating) && p.box.rating > 0)
    error('ch4_box_rule:rating', ...
          'p.box.rule is ''rating'' but p.box.rating is not a positive scalar.');
end

end
