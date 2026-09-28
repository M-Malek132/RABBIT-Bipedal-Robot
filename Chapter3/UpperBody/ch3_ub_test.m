function ch3_ub_test(gait_file)
%CH3_UB_TEST  Checks for the upper-body wrench extension.
%
%   ch3_ub_test()                          uses Results/ch3_baseline_n41.mat
%   ch3_ub_test('Results/<gait>.mat')
%
%   1  p.ub absent / 'none' reproduces f, g, lambda BIT FOR BIT.
%   2  J_sh (complex step) matches a central finite difference of Tt.
%   3  Q_ub(4:7) == 0: the wrench acts on x, z, q_t only.
%   4  Force through the hip (l_sh = 0) makes no torso moment: Q(3) == 0.
%   5  Virtual power: Q' dq == F . v_sh + M dqt.
%   6  Newton balance: raising the torso's vertical force by dF (torques
%      frozen) changes m_tot*a_com_z - GRF_z by exactly dF.
%   7  Pack/unpack/effective_params round-trip with optimize = true.

if nargin < 1, gait_file = 'Results/ch3_baseline_n41.mat'; end
S = load(gait_file);
p0 = S.p;
[X, ~, alpha] = ch3_col_unpack(S.z, p0);
x  = X(:, round(size(X, 2)/2));
q  = x(1:7);  dq = x(8:14);
tol = 1e-9;

% 1 ----------------------------------------------------------------------
[f0, g0, a0] = ch3_control_affine(x, p0);
p1 = ch3_ub_defaults(p0, 'none');
[f1, g1, a1] = ch3_control_affine(x, p1);
assert(isequal(f0, f1) && isequal(g0, g1) && isequal(a0.lam_drift, a1.lam_drift), ...
       'T1: mode none changed the dynamics');
fprintf('T1 pass: no wrench -> dynamics identical\n');

% 2 ----------------------------------------------------------------------
p2 = ch3_ub_defaults(p0, 'profile', 'fun', @(s) [30; -20; 5]);
[~, ~, ~, J] = ch3_ub_wrench(q, 0.5, p2);
Jfd = zeros(2, 7); h = 1e-6; pt = [0; -p2.ub.l_sh; 0; 1];
for k = 1:7
    e = zeros(7,1); e(k) = h;
    rp = Tt(q+e)*pt; rm = Tt(q-e)*pt;
    Jfd(:,k) = ([rp(1); rp(3)] - [rm(1); rm(3)]) / (2*h);
end
assert(max(abs(J(:) - Jfd(:))) < 1e-7, 'T2: J_sh mismatch');
fprintf('T2 pass: J_sh complex step vs FD, max err %.2e\n', max(abs(J(:)-Jfd(:))));

% 3 ----------------------------------------------------------------------
Q = ch3_ub_wrench(q, 0.5, p2);
assert(all(Q(4:7) == 0), 'T3: wrench leaked into leg coordinates');
fprintf('T3 pass: Q_ub(4:7) = 0, Q_ub(1:3) = [%.2f %.2f %.2f]\n', Q(1:3));

% 4 ----------------------------------------------------------------------
p4 = ch3_ub_defaults(p0, 'profile', 'fun', @(s) [30; -20; 0], 'l_sh', 0);
Q4 = ch3_ub_wrench(q, 0.5, p4);
assert(abs(Q4(3)) < tol, 'T4: force through hip made a moment');
fprintf('T4 pass: force at hip -> Q(3) = %.1e\n', Q4(3));

% 5 ----------------------------------------------------------------------
[Q5, W5, ~, J5] = ch3_ub_wrench(q, 0.5, p2);
lhs = Q5.' * dq;  rhs = W5(1:2).' * (J5*dq) + W5(3)*dq(3);
assert(abs(lhs - rhs) < tol, 'T5: virtual power mismatch');
fprintf('T5 pass: virtual power %.6f == %.6f\n', lhs, rhs);

% 6 ----------------------------------------------------------------------
% Total vertical momentum rate = sum_i m_i a_iz = lambda_z + Fz - m_tot g.
% With M ddq: the generalized-momentum identity along the pz direction
% gives, for u = 0,   m_tot * a_com_z = lambda_z + F_z - m_tot*g.
% We check the EQUIVALENT linear statement between two wrenches: raising Fz
% by dF (u fixed) changes  m_tot*a_com_z - lambda_z  by exactly dF.
m_tot = 74;
e2 = zeros(7,1); e2(2) = 1;   % pz direction: rigid vertical translation
pa = ch3_ub_defaults(p0, 'profile', 'fun', @(s) [0; 0; 0]);
pb = ch3_ub_defaults(p0, 'profile', 'fun', @(s) [0; 100; 0]);
[fa, ~, aa] = ch3_control_affine(x, pa);
[fb, ~, ab] = ch3_control_affine(x, pb);
% rigid-translation momentum: e2' M ddq = sum m_i * (d r_i/dpz) . a_i = -m_tot*a_com_z
pa_rate = e2.' * aa.M * fa(8:14);
pb_rate = e2.' * ab.M * fb(8:14);
d_mom = -(pb_rate - pa_rate);          % change of m_tot * a_com_z (pz is down)
d_lam = ab.lam_drift(2) - aa.lam_drift(2);
fprintf('T6: dFz = 100 N -> d(m a_com_z) = %.4f N, d(GRF_z) = %.4f N, sum check %.2e\n', ...
        d_mom, d_lam, abs(d_mom - d_lam - 100));
assert(abs(d_mom - d_lam - 100) < 1e-6, 'T6: vertical Newton balance violated');
fprintf('T6 pass: extra torso force = extra momentum rate - extra GRF\n');

% 7 ----------------------------------------------------------------------
B  = [10 20 30 20 10; 50 60 70 60 50; 1 2 3 2 1];
p7 = ch3_ub_defaults(p0, 'bezier', 'beta', B, 'optimize', true);
[X7, T7, al7] = ch3_col_unpack(S.z, p0);
if isfield(p0,'free_theta') && ~isempty(p0.free_theta) && p0.free_theta
    z7 = ch3_col_pack(X7, T7, al7, p7, [p0.theta_minus; p0.theta_plus]);
else
    z7 = ch3_col_pack(X7, T7, al7, p7);
end
assert(numel(z7) == numel(S.z) + numel(B), 'T7: pack length');
[X8, T8, al8] = ch3_col_unpack(z7, p7);
assert(isequal(X8, X7) && T8 == T7 && isequal(al8, al7), 'T7: unpack mismatch');
p8 = p7; p8.ub.beta = zeros(size(B));
p8 = ch3_col_effective_params(z7, p8);
assert(isequal(p8.ub.beta, B), 'T7: beta not decoded');
[lb, ub] = ch3_col_bounds(p7, size(X7, 2));
assert(numel(lb) == numel(z7) && numel(ub) == numel(z7), 'T7: bounds length');
fprintf('T7 pass: pack/unpack/effective_params/bounds consistent (+%d vars)\n', numel(B));

fprintf('\nAll upper-body tests passed.\n');
end
