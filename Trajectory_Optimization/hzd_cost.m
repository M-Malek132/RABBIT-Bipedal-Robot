%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST — penalizes only genuine simulation failure (impact event never
% fired, or NaNs). Equality/inequality constraints are handled by
% fmincon's nonlcon (hzd_constraints), not duplicated here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [~, total_torque_sq, ~, status] = simulate_hzd_gait(coeffs, x_start, p);

    if status < 0
        J = 1e6;   % failed / incomplete step — steer optimizer away
    else
        J = total_torque_sq;
    end
end
