function [ddq, tau, lambda] = hzd_control_and_dynamics(q, dq, coeffs, theta_minus, theta_plus, p, foot_ref)
%HZD_CONTROL_AND_DYNAMICS  The control law + constrained dynamics for one state.
%
%   [ddq, tau, lambda] = hzd_control_and_dynamics(q, dq, coeffs, ...
%                                    theta_minus, theta_plus, p, foot_ref)
%
% Computes the HZD virtual-constraint PD torque tau, the resulting joint
% accelerations ddq, and the stance-foot contact force lambda = [F_x; F_z]
% (world frame, up-positive; the force the ground applies to the robot).
%
% This is exactly the guts of hzd_closed_loop_ode, factored out so that the
% torque and contact force can be recomputed from a trajectory AFTER a sim
% (see gait_forces.m) without duplicating the controller. hzd_closed_loop_ode
% calls this and only appends the running torque-cost integrand.

    if nargin < 7, foot_ref = []; end

    % NOTE: this function no longer clamps theta_dot / e / de / tau / ddq to
    % "debug safety net" bounds. Those clamps silently saturated the outputs at
    % ad-hoc magnitudes (+/-100, +/-10, +/-1000), which never bind on a healthy
    % gait but could hand the ODE/optimizer a physically WRONG-but-finite ddq on
    % a divergent evaluation, hiding the failure. Divergent states are now left
    % to blow up naturally; simulate_hzd_gait flags any non-finite result as a
    % failed step (status < 0), so failures are visible instead of masked.

    theta     = theta_of_q(q);
    s         = max(0, min(1, (theta - theta_minus) / (theta_plus - theta_minus)));
    ds_dtheta = 1/(theta_plus - theta_minus);
    dtheta_dt = dtheta_dq_of(q) * dq;                % theta_dot = c*dq

    act_idx = 4:7;
    y  = q(act_idx);
    dy = dq(act_idx);

    yd  = zeros(p.nu,1);
    dyd = zeros(p.nu,1);
    for i = 1:p.nu
        [b, db] = bspline_eval(coeffs(i,:), s, p);
        yd(i)  = b;
        dyd(i) = db * ds_dtheta * dtheta_dt;         % d(yd)/dt via chain rule
    end

    e  = y  - yd;    % virtual-constraint output error
    de = dy - dyd;   % and its time derivative

    tau = -p.Kp.*e - p.Kd.*de;   % fixed-gain PD virtual-constraint law

    % Baumgarte only if a pin is given AND gains exist and are nonzero.
    use_baumgarte = ~isempty(foot_ref) && isfield(p,'baumgarte_beta') && ...
                    (p.baumgarte_beta ~= 0 || p.baumgarte_alpha ~= 0);
    if use_baumgarte
        [ddq, lambda] = rabbit_constrained_dynamics(q, dq, tau, foot_ref, ...
                                          p.baumgarte_alpha, p.baumgarte_beta);
    else
        [ddq, lambda] = rabbit_constrained_dynamics(q, dq, tau);
    end
end
