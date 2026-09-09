function tag = rabbit_generate_case_dynamics(s, out_dir)
%RABBIT_GENERATE_CASE_DYNAMICS  Rederive M, V, G from scratch for a scaled robot.
%
%   tag = rabbit_generate_case_dynamics(s)
%   tag = rabbit_generate_case_dynamics(s, out_dir)
%
% Chapter 4's mass_scale perturbation used to be applied by multiplying the
% NOMINAL M(q), V(q,dq), G(q) by a scalar s inside ch4_control_affine. This
% function instead reruns the SAME symbolic Lagrangian derivation used to
% build the nominal M.m/V.m/G.m in the first place --
% rabbit_energy_model_generalized_Lagrange.m -- but with the scaled link
% masses/inertias of Mass_Properties_scaled(s) from the start, and emits a
% dedicated, independently-derived set of dynamics functions for that case:
%
%   M_<tag>.m   mass matrix
%   V_<tag>.m   Coriolis/centrifugal vector (dM/dt * dq - d/dq of kinetic energy)
%   G_<tag>.m   gravity vector
%
% with tag = sprintf('scale%03.0f', round(s*100)), e.g. s = 0.7 -> 'scale070'.
%
% WHAT IS AND ISN'T REDERIVED. The transform chain (Tt, T1..T4) and the foot
% positions (P_st, P_sw) are pure kinematics -- no mass appears in a rotation
% or translation -- so they are geometry only and are reused unchanged from
% this folder. Only the three mass-dependent quantities are rebuilt.
%
% WHY THIS SHOULD MATCH THE OLD s*M(q) SHORTCUT, AND WHY IT IS STILL WORTH
% DOING. Every entry of Mass_Properties_scaled(s) is s times the nominal
% value (mass, COM-frame rotational inertia, and the parallel-axis
% correction m*(...) alike, since both are linear in mass for a fixed
% geometry). M is a sum of traces that are linear in each link's 4x4
% pseudo-inertia, so M_case = s*M_nominal follows algebraically. G is linear
% in mass by the same argument. This function does not take that argument on
% faith -- it reruns the derivation independently for each case and
% ch4_test_model compares the result numerically against the analytic
% s*M(q)/s*V(...)/s*G(q) form.
%
% Inputs
%   s       : mass_scale, positive finite scalar (1 reproduces the nominal
%             M.m/V.m/G.m and is not regenerated -- ch4_control_affine
%             already defers s==1 to ch3_control_affine)
%   out_dir : where to write the case files (default: this file's folder)
%
% Output
%   tag : the case tag used in the generated filenames
%
% See also MASS_PROPERTIES_SCALED, RABBIT_ENERGY_MODEL_GENERALIZED_LAGRANGE,
%          CH4_CONTROL_AFFINE.

if nargin < 2 || isempty(out_dir)
    out_dir = fileparts(mfilename('fullpath'));
end
if ~(isscalar(s) && isfinite(s) && s > 0)
    error('rabbit_generate_case_dynamics:s', ...
          'mass_scale must be a positive finite scalar (got %s).', mat2str(s));
end

tag = sprintf('scale%03.0f', round(s*100));
fprintf('--- rabbit_generate_case_dynamics: s = %.4f  (tag = %s) ---\n', s, tag);

[m, com, I] = Mass_Properties_scaled(s);

%% --- 1. mass matrix, same trace construction as the nominal derivation ---
syms theta [7 1] real

Tt_th = Tt(theta);
T1_th = T1(theta);
T2_th = T2(theta);
T3_th = T3(theta);
T4_th = T4(theta);

Msym = sym('M', [7 7]);
for j = 1:7
    for k = 1:7
        At = trace(diff(Tt_th, theta(j)) * I(:,:,1) * diff(transpose(Tt_th), theta(k)));
        A1 = trace(diff(T1_th, theta(j)) * I(:,:,2) * diff(transpose(T1_th), theta(k)));
        A2 = trace(diff(T2_th, theta(j)) * I(:,:,3) * diff(transpose(T2_th), theta(k)));
        A3 = trace(diff(T3_th, theta(j)) * I(:,:,4) * diff(transpose(T3_th), theta(k)));
        A4 = trace(diff(T4_th, theta(j)) * I(:,:,5) * diff(transpose(T4_th), theta(k)));
        Msym(j,k) = At + A1 + A2 + A3 + A4;
    end
end
Msym = simplify(Msym);

M_file = fullfile(out_dir, ['M_' tag]);
clear_generated(['M_' tag]);
matlabFunction(Msym, 'File', M_file, 'Vars', {theta});
fprintf('  wrote %s.m\n', M_file);

%% --- 2. Coriolis vector, via the Christoffel symbols of Msym directly ----
% The nominal derivation gets here through a d/dt(M(a(t))) chain-rule trick
% (treat theta as a function of time, differentiate, substitute the function
% and its derivative back to q, Dq). That trick is fragile under a SINGLE
% simultaneous subs of both a(t) and diff(a(t),t): substituting the function
% to a t-INDEPENDENT value before the derivative substitution is applied
% collapses every diff(a_i(t),t) term to zero as a side effect (verified by
% hand here -- it silently zeroed three of the seven rows of V for this
% model). The standard closed form for the SAME quantity sidesteps the
% chain-rule detour entirely:
%
%       v_i = sum_jk  0.5*(dM_ij/dq_k + dM_ik/dq_j - dM_jk/dq_i) * dq_j*dq_k
%
% which is exactly dM/dt*dq - d(1/2 dq'M dq)/dq written out in Christoffel
% symbols, so it is the same physical quantity -- just computed without the
% intermediate function-of-time object that caused the substitution bug.
syms omega [7 1] real

Vsym = sym(zeros(7,1));
for i = 1:7
    vi = sym(0);
    for j = 1:7
        for k = 1:7
            Cijk = 0.5 * ( diff(Msym(i,j), theta(k)) + diff(Msym(i,k), theta(j)) ...
                          - diff(Msym(j,k), theta(i)) );
            if Cijk == 0, continue; end
            vi = vi + Cijk * omega(j) * omega(k);
        end
    end
    Vsym(i) = simplify(vi);
end

in = [theta; omega];

V_file = fullfile(out_dir, ['V_' tag]);
clear_generated(['V_' tag]);
matlabFunction(Vsym, 'File', V_file, 'Vars', {in});
fprintf('  wrote %s.m\n', V_file);

%% --- 3. gravity vector, from the potential energy -------------------------
g0vec = [0 0 -9.8062 0]';

Pt = m(1) * g0vec.' * Tt_th * com(:,1);
P1 = m(2) * g0vec.' * T1_th * com(:,2);
P2 = m(3) * g0vec.' * T2_th * com(:,3);
P3 = m(4) * g0vec.' * T3_th * com(:,4);
P4 = m(5) * g0vec.' * T4_th * com(:,5);

P_energy = -(Pt + P1 + P2 + P3 + P4);

Gsym = [diff(P_energy,theta1); diff(P_energy,theta2); diff(P_energy,theta3); ...
        diff(P_energy,theta4); diff(P_energy,theta5); diff(P_energy,theta6); ...
        diff(P_energy,theta7)];
Gsym = simplify(Gsym);

G_file = fullfile(out_dir, ['G_' tag]);
clear_generated(['G_' tag]);
matlabFunction(Gsym, 'File', G_file, 'Vars', {theta});
fprintf('  wrote %s.m\n', G_file);

fprintf('--- done: %s ---\n\n', tag);

end

% ---------------------------------------------------------------------------
function clear_generated(name)
%CLEAR_GENERATED  Drop any in-memory copy of a previously generated case
% function before matlabFunction overwrites its file, so a stale cached
% version is never picked up by a later str2func/direct call in this or a
% subsequent MATLAB session.
clear(name);
rehash path;
end
