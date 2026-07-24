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

    % NOTE on the removed "debug safety net" clamps: the theta_dot/e/de/tau
    % clamps (+/-100, +/-10, +/-100, +/-1000) were removed -- measured over the
    % reference gait their peaks are 2.4 / 1.2 / 12.2 / 238, so they never bind
    % and removing them is a numerical no-op that stops a divergent evaluation
    % from being silently saturated into a wrong-but-finite value.
    %
    % The ddq clamp below is DIFFERENT and is KEPT: peak |ddq| on the reference
    % gait is ~2657 rad/s^2 (swing hip/knee at step start, s~0), so it DOES bind
    % over roughly the first ~13% of the step. This is real, well-conditioned
    % dynamics (cond(KKT)~6e2), not a singularity -- the reference gait genuinely
    % commands very large swing accelerations. The saved gait library was
    % optimized against these clamped dynamics (removing the clamp breaks its
    % periodicity: 7.5e-3 -> 0.35), so the clamp is load-bearing. It should be
    % retired by re-optimizing the gaits under a torque/GRF limit that bounds the
    % accelerations physically; until then, leaving it keeps the library valid.
    % See Results/ warm-start seed and hzd_problem_setup p.enforce_torque.

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

    % LOAD-BEARING clamp (see NOTE at top), toggled by p.clamp_ddq so a
    % re-optimization can search under the TRUE (unclamped) dynamics. Absent or
    % true => clamp ON, so saved p-structs and the existing gait library keep
    % their behaviour; set p.clamp_ddq = false to optimize a clamp-free gait.
    if ~isfield(p,'clamp_ddq') || p.clamp_ddq
        ddq = max(-1000, min(1000, ddq));
    end
end
