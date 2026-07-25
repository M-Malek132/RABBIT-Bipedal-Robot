function C = ch3_compare_controllers(z, p, u_box, n_steps)
%CH3_COMPARE_CONTROLLERS  The stage-8 payoff, measured.
%
%   C = ch3_compare_controllers(z, p)
%   C = ch3_compare_controllers(z, p, u_box, n_steps)
%
% Runs the SAME optimized gait under each control law and reports what each
% one actually does to the torque and to the outputs.  This is the experiment
% that justifies the last two stages of the chapter, so it is worth being
% precise about what it shows.
%
% The interesting column is what happens under a torque limit TIGHTER than the
% gait needs:
%
%   'iolin_pd'   has no notion of a limit. It commands whatever the inversion
%                asks for. On hardware the actuator saturates, the realized
%                torque is not the commanded torque, and every stability
%                guarantee that assumed ydd = mu is void -- silently.
%   'clfqp'      asks for the SMALLEST mu that still certifies the convergence
%                rate, so it is usually far cheaper than PD, but it still has
%                no hard bound: nothing stops the required mu from being large.
%   'clfqp_con'  carries the limit as a constraint. The torque is inside the
%                box BY CONSTRUCTION, and when the box and the CLF condition
%                conflict, the slack delta relaxes the CLF rather than letting
%                the QP go infeasible.
%
% READ delta.  A nonzero delta is the controller telling you it could not meet
% the convergence rate within the actuator limit. Per Remark 3.2, the
% exponential guarantee holds while delta = 0 and only under extra sufficient
% conditions once it is not -- so delta is not a diagnostic to skim past, it is
% the boundary of the theory.
%
% Inputs
%   z       : collocation decision vector (an optimized gait)
%   p       : parameter struct
%   u_box   : torque limit to impose for the comparison. Default: 60% of the
%             gait's own peak feedforward torque, so the box genuinely binds.
%   n_steps : steps to simulate per controller (default 4)
%
% Output
%   C : 1xk struct array, one entry per controller, with
%         .name .steps_completed .peak_torque .box_violation
%         .max_eta .delta_max .frac_delta_active .qp_infeasible .reason
%
% See also CH3_CTRL_CLF_QP, CH3_REPORT.

if nargin < 4 || isempty(n_steps), n_steps = 4; end

E = ch3_col_eval(z, p);
[X, ~, alpha] = ch3_col_unpack(z, p);
x0 = X(:,1);

peak_ff = max([abs(E.u(:)); abs(E.um(:))]);
if nargin < 3 || isempty(u_box)
    u_box = 0.6 * peak_ff;
end

names = {'iolin_pd', 'clfqp', 'clfqp_con'};
C = struct('name', {}, 'steps_completed', {}, 'peak_torque', {}, ...
           'box_violation', {}, 'max_eta', {}, 'delta_max', {}, ...
           'frac_delta_active', {}, 'qp_infeasible', {}, 'reason', {});

fprintf('\n%s\n CONTROLLER COMPARISON on one gait\n', repmat('=',1,74));
fprintf(' gait peak feedforward torque = %.1f Nm;  imposed box = %.1f Nm\n', ...
        peak_ff, u_box);
fprintf('%s\n', repmat('-',1,74));
fprintf(' %-11s %6s %10s %10s %10s %10s %8s\n', ...
        'controller', 'steps', 'peak |u|', 'over box', 'max|eta|', 'max delta', 'QP fail');
fprintf('%s\n', repmat('-',1,74));

for i = 1:numel(names)
    pc = p;
    pc.controller = names{i};
    pc.limits.u_max = u_box;
    % Only the constrained QP is TOLD about the box -- that is the point.
    pc.limits.enable.torque = strcmp(names{i}, 'clfqp_con');

    sim = ch3_simulate(x0, alpha, pc, n_steps);

    if sim.n_ok > 0
        F = ch3_forces(sim.t, sim.x, alpha, pc);
        peak   = F.torque_max;
        maxeta = max(vecnorm([F.y; F.ydot], 2, 1));
        dmax   = F.delta_max;
        fdel   = mean(F.delta > 1e-8);
        nqp    = F.qp_infeasible;
    else
        peak = NaN; maxeta = NaN; dmax = NaN; fdel = NaN; nqp = NaN;
    end

    C(i) = struct('name', names{i}, 'steps_completed', sim.n_ok, ...
                  'peak_torque', peak, 'box_violation', max(peak - u_box, 0), ...
                  'max_eta', maxeta, 'delta_max', dmax, ...
                  'frac_delta_active', fdel, 'qp_infeasible', nqp, ...
                  'reason', sim.reason);

    fprintf(' %-11s %6d %10.1f %10.1f %10.3e %10.2e %8d\n', ...
            names{i}, sim.n_ok, peak, C(i).box_violation, maxeta, dmax, nqp);
end

fprintf('%s\n', repmat('-',1,74));
fprintf(' "over box" is how far the COMMANDED torque exceeds the limit. Only\n');
fprintf(' clfqp_con is constrained to keep it at zero; the other two are not\n');
fprintf(' told the limit exists, which is precisely the gap stage 8 closes.\n');
fprintf('%s\n\n', repmat('=',1,74));

end
