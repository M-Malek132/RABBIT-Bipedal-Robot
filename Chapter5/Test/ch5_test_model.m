function pass = ch5_test_model()
%CH5_TEST_MODEL  The two plants, and the relative degrees the chapter rests on.
%
% Nothing above the model layer can be right if the relative degree is wrong,
% and a wrong relative degree is quiet: the barrier would still be built, the
% QP would still solve, and the enforced condition would simply be on the wrong
% derivative. So the relative degrees are MEASURED here, from the model, rather
% than read out of the parameter file.
%
% Checks:
%   1. springmass A, B reproduce eq (5.32)-(5.34) directly
%   2. the springmass is relative degree 6 from u to x3, and 1 to xdot1
%   3. L_f^6 y and L_g L_f^5 y against a hand derivation at unit parameters
%   4. the pendulum x0 of Fig. 5.4 is a genuine equilibrium
%   5. the pendulum is relative degree 4, and the decoupling matrix is D^-1 k/Jm
%   6. energy conservation with k -> stiff, xi = 0 and no input
%   7. eta is a change of coordinates: r*ny == nx on both plants
%
% See also CH5_SPRINGMASS, CH5_PENDULUM, CH5_TEST_BARRIER.

fprintf('\n=== ch5_test_model ===\n');
pass = true;

%% ------------------------------------------------------- 1. springmass A, B
p  = ch5_params();
sm = ch5_springmass(p);
m  = p.plant.m; k = p.plant.k;

% Written out from (5.32)-(5.34) independently of ch5_springmass's assembly.
A_ref = zeros(6); A_ref(1:3,4:6) = eye(3);
A_ref(4,1) = -k/m(1);  A_ref(4,2) =  k/m(1);
A_ref(5,1) =  k/m(2);  A_ref(5,2) = -2*k/m(2); A_ref(5,3) = k/m(2);
A_ref(6,2) =  k/m(3);  A_ref(6,3) = -k/m(3);
B_ref = [0;0;0; 1/m(1); 0; 0];

pass = rep('A matches (5.32)-(5.34)', norm(sm.A - A_ref, inf), 0, pass);
pass = rep('B matches (5.35)',        norm(sm.B - B_ref, inf), 0, pass);

% A uniform translation relaxes every spring, so it must be in null(A).
pass = rep('null(A) contains [1 1 1 0 0 0]', ...
           norm(sm.A * [1;1;1;0;0;0], inf), 1e-14, pass);

%% ------------------------------------------------ 2. springmass rel. degrees
b6 = ch5_barrier(ch5_x0(p), p);
ok = (b6.rb == 6);
fprintf('  [%s] %-34s rb = %d (want 6), L_g h = %g\n', tf(ok), ...
        'x3 <= x3max is rel. degree 6', b6.rb, b6.Lgh);
pass = pass && ok;

pv1 = ch5_params('constraint', struct('type','v1_max','value',1), ...
                 'ecbf.poles', 2.0);
b1 = ch5_barrier(ch5_x0(pv1), pv1);
ok = (b1.rb == 1) && abs(b1.Lgh + 1/m(1)) < 1e-14;
fprintf('  [%s] %-34s rb = %d (want 1), L_g h = %g\n', tf(ok), ...
        'xdot1 <= v1max is rel. degree 1', b1.rb, b1.Lgh);
pass = pass && ok;

%% ---------------------------------------------- 3. hand-derived Lie terms
% Expanding (5.32)-(5.34) by hand at k = m_i = 1:
%       x3^(6) = u - 4 x1 + 9 x2 - 5 x3
% so L_g L_f^5 y = 1 and L_f^6 y = -4x1 + 9x2 - 5x3. Derived on paper, checked
% here -- this is the one place the linear algebra of ch5_lie_rows is compared
% against something that did not come out of a matrix power.
if isequal(m, [1 1 1]) && k == 1
    rng(4);
    e_lf = 0; e_lg = 0;
    for j = 1:6
        xt = randn(6,1);
        io = ch5_io_lin(xt, p);
        e_lf = max(e_lf, abs(io.Lfry - (-4*xt(1) + 9*xt(2) - 5*xt(3))));
        e_lg = max(e_lg, abs(io.LgLfr1y - 1));
    end
    pass = rep('L_f^6 y == hand derivation', e_lf, 1e-12, pass);
    pass = rep('L_g L_f^5 y == k^2/(m1m2m3)', e_lg, 1e-14, pass);
else
    fprintf('  [skip] hand derivation assumes unit m, k\n');
end

%% --------------------------------------------- 4. the pendulum x0 is at rest
q  = ch5_params('system','pendulum');
z0 = ch5_x0(q);
[f0, ~] = ch5_control_affine(z0, q);
pass = rep('pendulum x0 is an equilibrium', norm(f0, inf), 1e-12, pass);

%% ------------------------------------- 5. pendulum rel. degree and decoupling
bq = ch5_barrier(z0, q);
ok = (bq.rb == 4) && all(abs(bq.Lgh) < 1e-14);
fprintf('  [%s] %-34s rb = %d (want 4), L_g h = %s\n', tf(ok), ...
        'py2 >= p2min is rel. degree 4', bq.rb, mat2str(bq.Lgh));
pass = pass && ok;

rng(5);
e_dec = 0; e_py = 0; min_eig = inf;
for j = 1:6
    xt = [randn(4,1); 0.5*randn(4,1)];
    d  = ch5_pendulum(xt, q.plant.pv);
    io = ch5_io_lin(xt, q);
    bt = ch5_barrier(xt, q);

    % L_g L_f^3 theta = D^-1 k / Jm, derived by hand from (5.36)-(5.37)
    ref = d.D \ (q.plant.k * eye(2)) / q.plant.Jm;
    e_dec = max(e_dec, norm(io.LgLfr1y - ref, inf) / norm(ref, inf));

    % and the constraint's version is that premultiplied by dpy/dtheta
    e_py = max(e_py, norm(bt.LgLfrb1 - d.Jpy * ref, inf) / max(1, norm(ref, inf)));

    min_eig = min(min_eig, min(eig(io.LgLfr1y)));
end
pass = rep('L_g L_f^3 y == D^-1 k / Jm', e_dec, 1e-12, pass);
pass = rep('L_g L_f^3 h == Jpy D^-1 k/Jm', e_py, 1e-12, pass);

ok = min_eig > 0;
fprintf('  [%s] %-34s min eig = %.3e\n', tf(ok), ...
        'decoupling matrix stays SPD', min_eig);
pass = pass && ok;

%% ------------------------------------------------- 6. energy, xi = 0, tau = 0
% With no damping and no input the pendulum-plus-motor system is conservative.
% This exercises D, C and Gv TOGETHER -- a sign error in the Coriolis term that
% no individual check would catch shows up here as drifting energy.
qe = ch5_params('system','pendulum');
qe.plant.xi = 0;
qe.plant.pv = ch5_pend_pv(qe);      % repack: pv is what the model reads

x = [ -2.2; 0.6; -2.0; 0.5; 0.3; -0.2; 0.1; 0.4 ];
E0 = pend_energy(x, qe);
dt = 1e-4;
for i = 1:20000
    k1 = ch5_ode_rhs(0, x,          zeros(2,1), qe);
    k2 = ch5_ode_rhs(0, x+dt/2*k1,  zeros(2,1), qe);
    k3 = ch5_ode_rhs(0, x+dt/2*k2,  zeros(2,1), qe);
    k4 = ch5_ode_rhs(0, x+dt*k3,    zeros(2,1), qe);
    x  = x + dt/6*(k1+2*k2+2*k3+k4);
end
E1 = pend_energy(x, qe);
pass = rep('energy conserved (xi=0, tau=0)', abs(E1-E0)/abs(E0), 1e-9, pass);

%% ------------------------------------------------ 7. no zero dynamics either way
for pp = {p, q}
    s = pp{1}.sys;
    ok = (s.r * s.ny == s.nx);
    fprintf('  [%s] %-34s r*ny = %d, nx = %d\n', tf(ok), ...
            sprintf('%s: eta spans the state', s.name), s.r*s.ny, s.nx);
    pass = pass && ok;
end

fprintf('--- ch5_test_model: %s ---\n', tf(pass));

end

% ---------------------------------------------------------------------------
function E = pend_energy(x, p)
%PEND_ENERGY  Total mechanical energy: links + motors + coupling springs.
d  = ch5_pendulum(x, p.plant.pv);
pvq = p.plant;
T_link  = 0.5 * d.dtheta.' * d.D * d.dtheta;
T_motor = 0.5 * pvq.Jm * (d.dthetam.' * d.dthetam);
U_grav  = -pvq.m(1)*pvq.g*pvq.lc(1)*cos(x(1)) ...
          -pvq.m(2)*pvq.g*(pvq.l(1)*cos(x(1)) + pvq.lc(2)*cos(x(1)+x(2)));
U_spring = 0.5 * pvq.k * sum((d.theta - d.thetam).^2);
E = T_link + T_motor + U_grav + U_spring;
end

function ok = rep(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-34s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
