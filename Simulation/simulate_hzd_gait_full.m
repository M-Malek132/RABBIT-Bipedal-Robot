%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAYBACK VERSION — full trajectory for animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_out, x_out] = simulate_hzd_gait_full(coeffs, x_start, p)
    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;

    xi0 = [x_start; 0];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.01, ...
                  'Events', @impact_event_wrapper);

    foot_ref = P_st(x_start(1:nq));   % pin the stance foot here (Baumgarte)

    [t_out, XI] = ode45(@(t,xi) hzd_ode_rhs(t,xi,coeffs,theta_minus,theta_plus,p,foot_ref), ...
                         [0 p.T_max], xi0, opts);

    x_out = XI(:, 1:2*nq);
end