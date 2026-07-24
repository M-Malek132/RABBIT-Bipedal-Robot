function [x_settled, history] = settle_initial_state(coeffs, x0, p, nIter)
%SETTLE_INITIAL_STATE Poincare-map (fixed-point) iteration for a periodic seed.
%
% Simulates one step, applies impact + reset, and feeds the result back in
% as the new initial state. If the closed-loop gait is even marginally
% stable this converges toward the periodic orbit, producing an x0 that is
% self-consistent with the step-to-step map -- especially for the VELOCITIES,
% which are essentially impossible to guess by hand and are what the
% periodicity constraints stall on.
%
% history(k) = ||x_next - x|| over the periodic components (idx 2:14).
% A decreasing history means the gait is settling onto an orbit; a flat or
% growing one means the closed-loop gait is not stable at these gains, and
% the optimizer will have a hard time regardless of seed.

    if nargin < 4, nIter = 15; end

    x = x0;
    x(1) = 0;                 % px is the translation gauge (bounds fix it at 0)
    history = zeros(nIter,1);

    fprintf('\n=== Settling initial state (Poincare iteration) ===\n');
    for k = 1:nIter
        [x_end, ~, ~, status] = simulate_hzd_gait(coeffs, x, p);

        if status < 0
            warning('settle_initial_state: step %d failed (status<0); returning last good state.', k);
            history = history(1:max(k-1,1));
            break;
        end

        x_next    = rabbit_reset_map(rabbit_impact_map(x_end));
        x_next(1) = 0;        % keep the translation gauge pinned

        history(k) = norm(x_next(2:14) - x(2:14));
        fprintf('  iter %2d: ||x_next - x|| = %.6f   theta_start = %+.4f\n', ...
                k, history(k), theta_of_q(x_next(1:p.nq)));

        x = x_next;
    end

    x_settled = x;

    fprintf('settled theta_start = %+.4f rad\n', theta_of_q(x_settled(1:p.nq)));
    r = check_ground_validity(x_settled(1:p.nq), x_settled(p.nq+1:end));
    if ~r.is_valid
        fprintf('NOTE: settled state has ground violations:\n');
        disp(r.violations);
    end
    fprintf('===================================================\n\n');
end
