function V = ch3_validity(sim, p)
%CH3_VALIDITY  How much of a run a real robot could have walked.
%
%   V = ch3_validity(sim, p)
%
% The simulation integrates the stance foot as a PIN: it neither lifts off nor
% slips, whatever the controller asks of the ground. A run can therefore
% "complete" steps whose contact force is pulling the foot into the floor, or
% whose friction demand no floor supplies, and a step count alone would score
% those as walking. This scores them against the contact they assume.
%
% A SAMPLE IS INVALID when the stance contact force at a recorded point (the
% steps' .lambda: at every solver point, under the torque actually held there,
% from the model the plant was integrated with) has
%
%       Fz <= 0                         the foot would lift off, or
%       |Fx| / Fz > p.limits.mu_s       it would slip (0.4, Table 3.1 -- the
%                                       coefficient the constrained QP's own
%                                       friction rows use)
%
% and A RUN IS VALID UP TO ITS FIRST INVALID SAMPLE: valid_steps is the number
% of completed steps before the step that contains it. One sample ends it, by
% design -- a real foot that lifts or slides once is no longer on the orbit the
% rest of the simulation assumes, so nothing after it is evidence of walking.
% The count and fraction of invalid samples are reported alongside, so a run
% ended by a single blip can be told from one that fights the ground throughout.
%
% Chapter 3's ch3_step and Chapter 4's ch4_step both record .lambda under
% sampled control (p.control_dt > 0). A continuous-control run, or one saved
% before 2026-09-18, carries none: V.available is then false and every count
% is NaN.
%
% Inputs
%   sim : ch3_simulate or ch4_simulate output
%   p   : parameter struct (uses p.limits.mu_s)
%
% Output
%   V : struct
%         .available     false when the run carries no contact forces
%         .valid_steps   completed steps before the first invalid sample
%         .first_step    step containing the first invalid sample (NaN: none)
%         .first_t       its time within that step [s]
%         .first_kind    'lift-off' | 'slip' | ''
%         .n_invalid     invalid samples over the completed steps
%         .frac_invalid  as a fraction of all samples
%         .Fz_min        smallest normal force [N]
%         .mu_max        largest |Fx|/Fz where Fz > 0
%         .mu_s          the coefficient scored against
%
% See also CH3_STEP, CH3_SIMULATE, CH3_COMPARE_CONTROLLERS, CH4_VALIDITY.

mu_s = p.limits.mu_s;
V = struct('available', false, 'valid_steps', NaN, 'first_step', NaN, ...
           'first_t', NaN, 'first_kind', '', 'n_invalid', NaN, ...
           'frac_invalid', NaN, 'Fz_min', NaN, 'mu_max', NaN, 'mu_s', mu_s);

steps = sim.steps;
if isempty(steps)
    V.available = true;
    V.valid_steps = 0; V.n_invalid = 0; V.frac_invalid = 0;
    return;
end
if ~isfield(steps, 'lambda') || any(arrayfun(@(s) isempty(s.lambda), steps))
    return;
end

V.available = true;
V.valid_steps = numel(steps);
n_all = 0; n_bad = 0;
Fz_min = inf; mu_max = 0;

for k = 1:numel(steps)
    lam = steps(k).lambda;
    Fx  = lam(1, :);
    Fz  = lam(2, :);
    ok  = isfinite(Fx) & isfinite(Fz);            % t = 0 of a failed period
    lift = ok & (Fz <= 0);
    mu   = abs(Fx) ./ Fz;
    slip = ok & (Fz > 0) & (mu > mu_s);
    bad  = lift | slip;

    n_all = n_all + nnz(ok);
    n_bad = n_bad + nnz(bad);
    if any(ok), Fz_min = min(Fz_min, min(Fz(ok))); end
    if any(ok & Fz > 0), mu_max = max(mu_max, max(mu(ok & Fz > 0))); end

    if isnan(V.first_step) && any(bad)
        j = find(bad, 1);
        V.first_step  = k;
        V.first_t     = steps(k).t(j);
        V.valid_steps = k - 1;
        if lift(j), V.first_kind = 'lift-off'; else, V.first_kind = 'slip'; end
    end
end

V.n_invalid    = n_bad;
V.frac_invalid = n_bad / max(n_all, 1);
V.Fz_min       = Fz_min;
V.mu_max       = mu_max;

end
