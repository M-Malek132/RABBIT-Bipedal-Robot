%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST — integral of squared torque PER UNIT STEP LENGTH.
%
%   J = (1/L_step) * int ||u||^2 dt
%
% This is Westervelt et al. eq. (6.43) / the thesis eq. (3.5). The division
% by step length is NOT cosmetic: without it, short steps are cheap, so the
% optimizer is actively rewarded for shrinking the stride -- which is exactly
% the pathology we measured (swing foot landing early at s~0.69, leaving
% theta_end persistently ~0.14 short of theta_plus). Normalizing removes that
% degenerate incentive and makes the cost "energy per distance travelled".
%
% Equality/inequality constraints are handled by fmincon's nonlcon
% (hzd_constraints), not duplicated here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [x_end, total_torque_sq, ~, status] = simulate_hzd_gait(coeffs, x_start, p);

    if status < 0
        J = 1e6;   % failed / incomplete step — steer optimizer away
        return;
    end

    % Step length = horizontal distance from the stance foot to the swing
    % foot at impact (Westervelt's p^h_2(q0^-)).
    q_end = x_end(1:p.nq);
    sw    = P_sw(q_end);
    st    = P_st(q_end);
    L_step = sw(1) - st(1);

    % Normalize by step length, but CLAMP the divisor instead of rejecting
    % short steps outright.
    %
    % An earlier version returned a constant 1e6 whenever L_step fell below
    % p.min_step_length. That made the cost a FLAT PLATEAU: at the initial
    % guess (whose first step is shorter than the threshold) the gradient was
    % identically zero, fmincon had no descent direction, and it bolted into a
    % far worse region trying to escape. Clamping keeps the cost finite and
    % CONTINUOUS everywhere, so a gradient always exists.
    %
    % A short or backward step is then discouraged smoothly by the quadratic
    % penalty below rather than by a cliff. Step length is really pinned by
    % the theta_end equality constraint; this term only stops the optimizer
    % from drifting to a degenerate stride while it gets there.
    L_eff = max(L_step, p.min_step_length);
    J = total_torque_sq / L_eff;

    shortfall = max(0, p.min_step_length - L_step);
    J = J + p.short_step_weight * shortfall^2;
end
