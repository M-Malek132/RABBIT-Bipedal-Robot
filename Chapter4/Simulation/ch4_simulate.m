function out = ch4_simulate(x0, alpha, p, n_steps, xi0)
%CH4_SIMULATE  Integrate n_steps of the uncertain hybrid system.
%
%   out = ch4_simulate(x0, alpha, p, n_steps)
%   out = ch4_simulate(x0, alpha, p, n_steps, xi0)
%
% Chains ch4_step, carrying BOTH the robot state and the controller state from
% one step to the next. The controller state is what makes this more than a
% loop: the L1 estimates accumulate across footstrikes, so step 5 starts with
% four steps' worth of learning and its behaviour is not a repeat of step 1.
% That accumulation is the mechanism behind the Section 4.2.4 result that the
% rate of convergence is unchanged across model perturbations.
%
% ------------------------------------------------- the randomized-load study
% If p.load_random_range is nonempty, the carried torso mass is REDRAWN AT
% EVERY STEP from that interval -- the experiment behind Fig. 4.11a, where the
% load varies randomly in 0-30 kg from step to step and the controller is never
% told. The draw happens between steps, so within a step the load is constant
% and the plant is a well-defined time-invariant system; the discontinuity sits
% at the footstrike where the dynamics already jump.
%
% This is the sharpest test of the adaptive scheme in the chapter: a robust
% controller must be sized for the worst case of the whole interval on every
% step, while the adaptive one re-converges to whatever the current step
% actually carries.
%
% Inputs
%   x0      : 14x1 start state
%   alpha   : ny x n_ctrl coefficients
%   p       : parameter struct
%   n_steps : number of steps to attempt
%   xi0     : initial controller state, or [] / omitted to initialize
%
% A run stops at the first step that does not reach the guard OR that reaches
% it by falling (see not_a_step below); n_ok counts real steps only, and a
% rejected step contributes nothing to .t / .x, exactly like one that timed out.
%
% Outputs
%   out : struct with
%           .steps   1 x k struct array of ch4_step outputs
%           .t .x    concatenated trajectory (t monotone across steps)
%           .xi .t_xi   concatenated controller state and its sample grid
%           .loads   carried mass used on each ATTEMPTED step [kg]
%           .n_ok .failed .reason .x_final .xi_final
%
% See also CH4_STEP, CH3_SIMULATE, CH4_FORCES.

if nargin < 5, xi0 = []; end

steps  = {};
t_all  = [];
x_all  = [];
xi_all = [];
txi_all = [];
loads  = [];

t_off  = 0;
x      = x0(:);
xi     = xi0;
failed = false;
reason = 'completed';

randomize = ~isempty(p.load_random_range);
if randomize
    % Reproducible load sequence. p.load_seed picks which one; a parameter
    % struct saved before that field existed ran seed 11.
    seed = 11;
    if isfield(p, 'load_seed') && ~isempty(p.load_seed), seed = p.load_seed; end
    rng(seed);
end

% A STRUCTURED uncertainty (ch4_uncertainty_set) brings three things that
% outlive a single step: one measurement-noise stream for the whole run (a
% handle, so its state carries across steps), the command pending in a
% one-sample actuation delay, and -- when the two legs were perturbed
% independently -- which PHYSICAL leg is the stance leg. The model always calls
% the stance leg links 2-3 and joints 1-2, and the legs relabel at every
% impact, so on even steps the drawn parameters of the two legs are swapped.
S_unc = ch4_structured_part(p.uncertainty);
noise_rs = [];
if isfield(p.uncertainty, 'noise') && isstruct(p.uncertainty.noise) ...
        && ((isfield(p.uncertainty.noise, 'q')  && p.uncertainty.noise.q  > 0) || ...
            (isfield(p.uncertainty.noise, 'dq') && p.uncertainty.noise.dq > 0))
    nseed = 0;
    if isfield(p.uncertainty.noise, 'seed'), nseed = p.uncertainty.noise.seed; end
    noise_rs = RandStream('mt19937ar', 'Seed', nseed);
end
swap_legs = ~isempty(S_unc) && isfield(p.uncertainty, 'asymmetric') ...
            && ~isempty(p.uncertainty.asymmetric) && p.uncertainty.asymmetric;
u_prev = [];

for k = 1:n_steps

    if ~all(isfinite(x))
        failed = true;
        reason = sprintf('non-finite state entering step %d', k);
        break;
    end

    pk = p;
    if randomize
        lo = p.load_random_range(1);
        hi = p.load_random_range(2);
        pk.uncertainty.load_mass = lo + (hi - lo) * rand();
    end
    loads(end+1) = pk.uncertainty.load_mass; %#ok<AGROW>
    if ~isempty(noise_rs), pk.uncertainty.noise.stream = noise_rs; end
    if ~isempty(S_unc),    pk.uncertainty.u_prev = u_prev;         end
    if swap_legs && mod(k, 2) == 0
        pk.uncertainty = swap_leg_parameters(pk.uncertainty);
    end

    s = ch4_step(x, xi, alpha, pk);
    if isfield(s, 'u_last'), u_prev = s.u_last; end

    if ~s.ok
        failed = true;
        if isfield(s, 'contact_invalid') && s.contact_invalid
            % p.stop_on_invalid: the step ended where the ground could not
            % have held the pinned foot; it is not counted as walked.
            reason = sprintf(['step %d lost contact validity: %s at t = %.3f s ' ...
                              '(Fz %.1f N, |Fx|/Fz %.2f)'], k, s.invalid.kind, ...
                             s.invalid.t, s.invalid.Fz, s.invalid.mu);
        else
            reason = sprintf('step %d never reached the guard (T = %.3f s)', k, s.T);
        end
        break;
    end

    why = not_a_step(s, pk);
    if ~isempty(why)
        failed = true;
        reason = sprintf('step %d is not a step: %s', k, why);
        break;
    end

    steps{end+1} = s;                        %#ok<AGROW>
    t_all = [t_all, t_off + s.t];            %#ok<AGROW>
    x_all = [x_all, s.x];                    %#ok<AGROW>
    if ~isempty(s.xi)
        xi_all  = [xi_all,  s.xi];           %#ok<AGROW>
        txi_all = [txi_all, t_off + s.t_xi]; %#ok<AGROW>
    end

    t_off = t_off + s.T;
    x     = s.x_next;
    xi    = s.xi_next;
end

if isempty(steps)
    % An EMPTY STRUCT ARRAY WITH THE RIGHT FIELDS, not struct([]).
    %
    % struct([]) has no fields at all, so the perfectly ordinary idiom
    % [sim.steps.T] -- used by ch4_report, ch4_compare_controllers and any
    % analysis script -- raises "Unrecognized field name" instead of returning
    % []. A completely failed run is a normal outcome in Chapter 4, not an
    % exceptional one, so it must flow through the same code path as a good
    % run rather than blowing up in the reporting layer.
    steps_arr = struct('t', {}, 'x', {}, 'x_end', {}, 'T', {}, 'ok', {}, ...
                       'sol', {}, 'xi', {}, 't_xi', {}, 'xi_end', {}, ...
                       'L_step', {}, 'x_next', {}, 'impulse', {}, ...
                       'xi_next', {}, 'u', {}, 't_u', {}, 'lambda', {}, ...
                       'contact_invalid', {}, 'invalid', {}, 'u_last', {});
else
    steps_arr = [steps{:}];
end

out = struct('steps', steps_arr, 't', t_all, 'x', x_all, ...
             'xi', xi_all, 't_xi', txi_all, 'loads', loads, ...
             'n_ok', numel(steps), 'failed', failed, 'reason', reason, ...
             'x_final', x, 'xi_final', xi);

end

% ---------------------------------------------------------------------------
function unc = swap_leg_parameters(unc)
%SWAP_LEG_PARAMETERS  The other physical leg is the stance leg on this step.
unc.links = unc.links([1 4 5 2 3]);
for f = {'J_ref', 'b_visc', 'tau_coulomb', 'tau_bias'}
    if isfield(unc, f{1}) && numel(unc.(f{1})) == 4
        v = unc.(f{1});
        unc.(f{1}) = v([3 4 1 2]);
    end
end
end

% ---------------------------------------------------------------------------
function why = not_a_step(s, p)
%NOT_A_STEP  Reject a guard crossing that is a fall, not a footstrike.
%
% The guard is only "swing foot reaches the ground", and a falling robot
% reaches the ground too. Counting those crossings as steps is how a run that
% hit max||eta|| = 290 and put its hip through the floor was tabulated as
% "3 steps completed". Measured on such runs, every spurious crossing shows at
% least one of four signatures that no real step shows -- perturbed and
% badly tracked steps included, down to 0.159 s and 0.186 m, with the torso
% never past 0.47 rad:
%
%   * T at the guard's first armed instant (0.050 s): the swing foot was
%     already on or below the ground when p.guard_min_time expired, i.e. it
%     never left it;
%   * the hip below half its nominal height (seen: 0.035 m, then -0.146 m);
%   * the torso pitched past 1 rad (57 deg), Chapter 3's default torso pitch
%     box -- seen: a 0.54 s "step" lunging 0.93 m with the torso
%     at 1.07 rad, followed by a 0.15 s one in which it rotated on to 4.12 rad
%     (a somersault) while every other signature still passed;
%   * a step length under p.step_len_min -- the floor Chapter 3 imposes to
%     rule out stepping in place -- or negative (seen: +0.051, -0.692, -1.138).
%
% Returns '' for a real step, otherwise the reason.
TORSO_MAX = 1.0;                       % [rad], ch3_params' default qt_range
why = '';

if s.T <= p.guard_min_time + max(p.control_dt, 1e-3)
    why = sprintf(['the guard fired as soon as it was armed (T = %.3f s): ' ...
                   'the swing foot never left the ground'], s.T);
    return;
end

hip_min = min(-s.x(2, :));             % pz is DOWN-positive
if hip_min < 0.5 * p.limits.hip_h
    why = sprintf('the hip fell to %.3f m', hip_min);
    return;
end

[qt_worst, j] = max(abs(s.x(3, :)));
if qt_worst > TORSO_MAX
    why = sprintf('the torso pitched to %+.2f rad', s.x(3, j));
    return;
end

if s.L_step < p.step_len_min
    why = sprintf('step length %.3f m, under the %.2f m floor', ...
                  s.L_step, p.step_len_min);
end
end
