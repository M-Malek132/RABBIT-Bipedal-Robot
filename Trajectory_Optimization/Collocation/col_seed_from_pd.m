function z0 = col_seed_from_pd(coeffs, x_start, p)
%COL_SEED_FROM_PD  Warm-start the collocation from a PD-controlled step.
%
%   z0 = col_seed_from_pd(coeffs, x_start, p)
%
% Simulates ONE step under the fixed-gain PD virtual-constraint law (the
% existing, fast shooting sim), resamples the trajectory to p.N nodes uniform
% in time, and recomputes the node torques and contact forces from that same
% controller. Because those (q,dq,u,lambda) satisfy M ddq + C + G = B u +
% J_st' lambda by construction, the Hermite-Simpson defects and the EOM start
% nearly satisfied -- a good seed is what makes the collocation NLP converge.
%
% The virtual constraints y = q(4:7) - yd(s) will start slightly nonzero (PD
% has tracking error); fmincon drives them to the manifold.

    N  = p.N;   nq = p.nq;   nu = p.nu;   nx = 2*nq;

    p_pd = p;   p_pd.controller = 'pd';     % force the PD law for the seed sim
    [t, x] = simulate_hzd_gait_full(coeffs, x_start, p_pd);   % x: M x 14

    T  = t(end) - t(1);
    tn = linspace(t(1), t(end), N);
    Xn = interp1(t, x, tn).';                % nx x N  (interp1 works column-wise)

    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;

    U   = zeros(nu, N);
    Lam = zeros(2,  N);
    for k = 1:N
        q  = Xn(1:nq,k);
        v  = Xn(nq+1:nx,k);
        [~, tau, lam] = hzd_control_and_dynamics(q, v, coeffs, ...
                                theta_minus, theta_plus, p_pd);
        U(:,k)   = tau;
        Lam(:,k) = lam;
    end

    z0 = col_pack(Xn, U, Lam, coeffs, T);
end
