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
% Outputs
%   out : struct with
%           .steps   1 x k struct array of ch4_step outputs
%           .t .x    concatenated trajectory (t monotone across steps)
%           .xi .t_xi   concatenated controller state and its sample grid
%           .loads   1 x k carried mass actually used on each step [kg]
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
    rng(11);                              % reproducible load sequence
end

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

    s = ch4_step(x, xi, alpha, pk);

    if ~s.ok
        failed = true;
        reason = sprintf('step %d never reached the guard (T = %.3f s)', k, s.T);
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
                       'xi_next', {});
else
    steps_arr = [steps{:}];
end

out = struct('steps', steps_arr, 't', t_all, 'x', x_all, ...
             'xi', xi_all, 't_xi', txi_all, 'loads', loads, ...
             'n_ok', numel(steps), 'failed', failed, 'reason', reason, ...
             'x_final', x, 'xi_final', xi);

end
