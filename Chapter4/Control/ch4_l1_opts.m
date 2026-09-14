function o = ch4_l1_opts(p)
%CH4_L1_OPTS  The L1 design choices in p.l1, resolved in one place.
%
%   o = ch4_l1_opts(p)
%
% The L1 law reads its structural options -- which predictor, the predictor's
% error-injection rate, alpha's own adaptation gain and regressor cap, and
% whether the constrained law bounds the torque it actually applies -- through
% here rather than from p.l1 directly, for one reason: a parameter struct saved
% before these options existed (an old ch4_result .mat) does not carry them,
% and re-analysing it must reproduce the controller that run actually used. So
% a missing field resolves to the Section 4.2 formulation: the thesis
% predictor, one gain for both estimates, no cap, and the box on mu1 alone.
%
% Output
%   o : struct .predictor ('thesis' | 'plant') .predictor_rate [rad/s]
%              .Gamma .Gamma_alpha
%              .phi_max (cap on alpha's regressor; Inf for none)
%              .constrain_applied
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

% The alpha regressor's cap, phi_max, from the fastest the estimator loop may
% run in radians per control sample (see ch4_l1_deriv). The loop's natural
% frequency is sqrt(Gamma + Gamma_alpha phi^2), so holding it at kappa/dt gives
% phi_max = sqrt((kappa/dt)^2 - Gamma) / sqrt(Gamma_alpha). No cap for
% continuous control, a frozen alpha, or kappa unset -- the law as written.
o.phi_max = Inf;
kappa = 0;
if isfield(p.l1, 'alpha_regressor_rate') && ~isempty(p.l1.alpha_regressor_rate)
    kappa = p.l1.alpha_regressor_rate;
end
if kappa > 0 && p.control_dt > 0 && o.Gamma_alpha > 0
    room = (kappa / p.control_dt)^2 - o.Gamma;
    if room <= 0
        error('ch4_l1_opts:regressorRate', ...
              ['p.l1.alpha_regressor_rate = %g rad per sample leaves no room for ' ...
               'alpha: the beta channel alone already runs at %.2f rad per ' ...
               'sample.'], kappa, sqrt(o.Gamma) * p.control_dt);
    end
    o.phi_max = sqrt(room / o.Gamma_alpha);
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
