function out = ch3_simulate(x0, alpha, p, n_steps)
%CH3_SIMULATE  Integrate n_steps of the hybrid system.
%
%   out = ch3_simulate(x0, alpha, p, n_steps)
%
% Chains ch3_step, seeding each step from the previous step's post-Delta
% state.  Stops early and reports if a step fails to reach the guard or the
% state stops being finite -- which is how a periodic-but-UNSTABLE gait
% announces itself: the velocities grow step over step until the integration
% blows up.  Periodicity alone does not imply walking; see ch3_poincare.
%
% Inputs
%   x0      : 14x1 start state
%   alpha   : ny x n_ctrl coefficients
%   p       : parameter struct
%   n_steps : number of steps to attempt
%
% Outputs
%   out : struct with
%           .steps   1 x k struct array of ch3_step outputs
%           .t       concatenated time (monotone across steps)
%           .x       concatenated state trajectory
%           .n_ok    number of completed steps
%           .failed  true if it stopped early
%           .reason  why it stopped
%
% See also CH3_STEP, CH3_POINCARE.

% Collected as a cell array, then concatenated. Pre-declaring the struct
% fields here instead would silently couple this file to ch3_step's exact
% field list, and adding a field there would break assignment here.
steps = {};

t_all = [];
x_all = [];
t_off = 0;
x     = x0(:);
failed = false;
reason = 'completed';

for k = 1:n_steps
    if ~all(isfinite(x))
        failed = true; reason = sprintf('non-finite state entering step %d', k);
        break;
    end

    s = ch3_step(x, alpha, p);

    if ~s.ok
        failed = true;
        reason = sprintf('step %d never reached the guard (T = %.3f s)', k, s.T);
        break;
    end

    steps{end+1} = s;               %#ok<AGROW>
    t_all = [t_all, t_off + s.t];   %#ok<AGROW>
    x_all = [x_all, s.x];           %#ok<AGROW>
    t_off = t_off + s.T;
    x = s.x_next;
end

if isempty(steps)
    steps_arr = struct([]);
else
    steps_arr = [steps{:}];
end

out = struct('steps', steps_arr, 't', t_all, 'x', x_all, ...
             'n_ok', numel(steps), 'failed', failed, 'reason', reason, ...
             'x_final', x);

end
