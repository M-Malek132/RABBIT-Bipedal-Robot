%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIT SIMULATION — now with a REAL status flag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_end, total_torque_sq, max_penetration, status] = simulate_hzd_gait(coeffs, x_start, p)
    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;

    xi0 = [x_start; 0];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.01, ...
                  'Events', @impact_event_wrapper);

    sol = ode45(@(t,xi) hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p), ...
                [0 p.T_max], xi0, opts);

    xi_end          = sol.y(:,end);
    x_end           = xi_end(1:2*nq);
    total_torque_sq = xi_end(end);

    % STATUS: valid only if the swing foot actually triggered the impact
    % event (as opposed to running out of T_max with the foot never
    % landing), and the result contains no NaNs.
    impacted = isfield(sol,'ie') && ~isempty(sol.ie);
    hit_nan  = any(isnan(xi_end));
    if impacted && ~hit_nan
        status = 1;
    else
        status = -1;
    end

    % Dense-sample the whole step to find worst ground clipping.
    % World-frame Z is UP-positive (ground = 0), so penetration depth is
    % -min(height); a negative max_penetration means everyone stayed
    % above ground by that margin.
    t_samples  = linspace(sol.x(1), sol.x(end), 40);
    xi_samples = deval(sol, t_samples);

    max_penetration = -inf;
    for k = 1:length(t_samples)
        q_k = xi_samples(1:nq, k);
        pts = get_body_points(q_k);
        h_min = min([pts.torso_top(2), pts.stance_knee(2), pts.swing_knee(2), pts.swing_foot(2)]);
        max_penetration = max(max_penetration, -h_min);
    end
end
