function [t_out, x_out, torque_cost] = simulate_clf_gait(coeffs, x_start, p)
%SIMULATE_CLF_GAIT  One single-support step under the CLF-QP controller.
%
%   [t_out, x_out, torque_cost] = simulate_clf_gait(coeffs, x_start, p)
%
% Convenience wrapper: forces p.controller = 'clfqp' and runs the ordinary
% single-support simulation, so the CLF-QP law goes through the SAME dispatched
% path (hzd_ode_rhs) and event/impact machinery as everything else -- there is
% no separate integrator to drift out of sync. Returns the state trajectory and
% the integrated torque-squared cost.
%
% NOTE: a QP is solved at every ODE evaluation, so this is slower than the PD
% sim. Meant for validating / comparing the CLF-QP law on a known gait.

    p.controller = 'clfqp';

    [t_out, x_out] = simulate_hzd_gait_full(coeffs, x_start, p);

    % Only pay for the (separate) scalar-cost sim if the caller wants the cost;
    % each sim solves a QP per ODE step, so avoid running it twice needlessly.
    if nargout > 2
        [~, torque_cost] = simulate_hzd_gait(coeffs, x_start, p);
    end
end
