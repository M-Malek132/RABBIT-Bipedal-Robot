function o = ch4_l1_opts(p)
%CH4_L1_OPTS  The L1 design choices in p.l1, resolved in one place.
%
%   o = ch4_l1_opts(p)
%
% The L1 law reads its structural options -- which predictor, the predictor's
% error-injection rate, alpha's own adaptation gain, and whether the
% constrained law bounds the torque it actually applies -- through here rather
% than from p.l1 directly, for one reason: a parameter struct saved before
% these options existed (an old ch4_result .mat) does not carry them, and
% re-analysing it must reproduce the controller that run actually used. So a
% missing field resolves to the Section 4.2 formulation: the thesis predictor,
% one gain for both estimates, and the box on mu1 alone.
%
% Output
%   o : struct .predictor ('thesis' | 'plant') .predictor_rate [rad/s]
%              .Gamma .Gamma_alpha .constrain_applied
%
% See also CH4_PARAMS, CH4_L1_DERIV, CH4_L1_ADVANCE, CH4_CTRL_L1.

o = struct('predictor', 'thesis', 'predictor_rate', 0, ...
           'Gamma', p.l1.Gamma, 'Gamma_alpha', p.l1.Gamma, ...
           'constrain_applied', false);

if isfield(p.l1, 'predictor') && ~isempty(p.l1.predictor)
    o.predictor = lower(p.l1.predictor);
end
if isfield(p.l1, 'predictor_rate') && ~isempty(p.l1.predictor_rate)
    o.predictor_rate = p.l1.predictor_rate;
end
if isfield(p.l1, 'Gamma_alpha') && ~isempty(p.l1.Gamma_alpha)
    o.Gamma_alpha = p.l1.Gamma_alpha;
end
if isfield(p.l1, 'constrain_applied') && ~isempty(p.l1.constrain_applied)
    o.constrain_applied = logical(p.l1.constrain_applied);
end

switch o.predictor
    case 'thesis'
    case 'plant'
        if ~(isfinite(o.predictor_rate) && o.predictor_rate > 0)
            error('ch4_l1_opts:rate', ...
                  ['p.l1.predictor_rate must be positive for the plant ' ...
                   'predictor (got %g): it is the rate at which the ' ...
                   'prediction error decays, and at zero the error never ' ...
                   'forgets its initial value.'], o.predictor_rate);
        end
    otherwise
        error('ch4_l1_opts:predictor', ...
              'Unknown p.l1.predictor "%s" (expected thesis|plant).', ...
              o.predictor);
end

end
