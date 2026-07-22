function F = gait_forces(coeffs, x_start, p)
%GAIT_FORCES  Physical quantities of a gait for the Table 3.1 constraints.
%
%   F = gait_forces(coeffs, x_start, p)
%
% DIAGNOSTIC ONLY -- nothing is enforced here. Re-simulates the step and
% recomputes, at every sample, the joint torque and stance-foot contact force
% (via hzd_control_and_dynamics), plus the impact impulse (via
% rabbit_impact_map). These are the raw numbers needed to choose physically
% sensible limits for THIS robot before adding them as constraints.
%
% RABBIT here is ~30 kg (weight ~294 N) with NO gear ratio (u is the JOINT
% torque). The thesis Table 3.1 values are ATRIAS-specific (63 kg, 50:1
% drives) and do not transfer directly -- use these measurements instead.
%
% Fields (all SI):
%   F.peak_torque         [nu x 1] max |u_i| over the step, per joint [Nm]
%   F.peak_torque_max     scalar   max over joints [Nm]
%   F.min_vertical_grf    scalar   min vertical stance GRF over the step [N]
%   F.mean_vertical_grf   scalar   mean vertical stance GRF [N] (~body weight)
%   F.max_friction_ratio  scalar   max |F_x / F_z| over the step (want < mu)
%   F.impact_impulse      scalar   |impact impulse| at touchdown [Ns]
%   F.impact_impulse_vec  [2 x 1]  the impulse [Ix; Iz] [Ns]
%
% NOTE ON SIGN: F_z is the vertical component of the contact multiplier. A
% healthy stance has the ground pushing UP, so the mean should be ~ +weight.
% If it comes out consistently negative, the KKT sign convention is flipped
% and the vertical GRF is -F_z (the ratio |F_x/F_z| is sign-agnostic).

    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;
    foot_ref    = P_st(x_start(1:nq));   % same pin the simulator uses

    [~, x] = simulate_hzd_gait_full(coeffs, x_start, p);
    N = size(x,1);

    peak_torque = zeros(p.nu,1);
    Fz_min  =  inf;
    Fz_sum  =  0;
    mu_max  =  0;

    for k = 1:N
        q  = x(k,1:nq)';
        dq = x(k,nq+1:end)';
        [~, tau, lambda] = hzd_control_and_dynamics(q, dq, coeffs, ...
                                theta_minus, theta_plus, p, foot_ref);
        peak_torque = max(peak_torque, abs(tau));

        Fx = lambda(1);   Fz = lambda(2);   % tangential, vertical (world)
        Fz_min = min(Fz_min, Fz);
        Fz_sum = Fz_sum + Fz;
        % Friction cone only makes physical sense where the foot is actually
        % in compression (Fz > 0). Where Fz crosses zero the ratio |Fx/Fz|
        % blows up -- that is a GROUND-REACTION problem (foot lifting), not a
        % friction one, and is caught by the min-GRF constraint instead. Only
        % sample the ratio where Fz exceeds a small positive threshold.
        if Fz > 1.0   % [N]
            mu_max = max(mu_max, abs(Fx) / Fz);
        end
    end

    % Impact impulse at touchdown (end of the step)
    [~, impulse] = rabbit_impact_map(x(end,:)');

    F.peak_torque        = peak_torque;
    F.peak_torque_max    = max(peak_torque);
    F.min_vertical_grf   = Fz_min;
    F.mean_vertical_grf  = Fz_sum / max(N,1);
    F.max_friction_ratio = mu_max;
    F.impact_impulse     = norm(impulse);
    F.impact_impulse_vec = impulse;
end
