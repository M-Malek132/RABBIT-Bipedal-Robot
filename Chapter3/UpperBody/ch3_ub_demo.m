%CH3_UB_DEMO  Effect of an upper-body wrench on RABBIT, in three parts.
%
%   Run from the repository root after startup.
%
%   PART A  Sensitivity on the nominal orbit (seconds, no optimizer).
%           Holds the gait FIXED and asks: what extra joint torque and ground
%           reaction does the wrench cost if the legs must keep walking exactly
%           as before?  u_ff is re-solved at every collocation node with the
%           wrench in the drift, so this is exact on the orbit.
%
%   PART B  Closed loop (the saved gait's controller, alpha unchanged).
%           Lets the wrench actually push the robot and shows the response of
%           x (progression), z (hip height) and q_t over several steps.
%
%   PART C  Re-optimize the gait WITH the wrench (needs the Optimization
%           Toolbox). Two sub-cases:
%             C1  fixed profile   : new alpha for a known upper-body load
%             C2  optimized wrench: beta solved with alpha, crutch-like
%                                   (Fz >= 0), penalized by p.ub.w
%
% The test profile used in A-C1 is ILLUSTRATIVE ONLY -- a one-cycle-per-step
% sagittal pattern of plausible magnitude. Replace ub_fun with measured data
% (e.g. shoulder/crutch reaction forces from an OpenSim model, resampled on
% the phase s) before drawing conclusions from it.
%
% See also CH3_UB_DEFAULTS, CH3_UB_WRENCH, CH3_UB_TEST.

gait_file = fullfile('Results', 'ch3_gait_posture_195.mat');   % Chapter 4's default gait
n_sim     = 8;                                                   % steps in part B
do_optim  = exist('fmincon', 'file') == 2;                       % part C needs fmincon

% Illustrative upper-body wrench, as a function of the phase s in [0,1]:
%   Fx : +-20 N fore-aft, Fz : +-30 N vertical, M : +-10 N m torso pitch
ub_fun = @(s) [20*sin(2*pi*s); 30*sin(4*pi*s); 10*cos(2*pi*s)];

S  = load(gait_file);
p0 = ch3_upgrade_params(S.p);
[X, T, alpha] = ch3_col_unpack(S.z, p0);
N  = size(X, 2);
pW = ch3_ub_defaults(p0, 'profile', 'fun', ub_fun);

%% ------------------------------------------------------------ PART A
u0 = zeros(p0.nu, N);  uW = u0;  L0 = zeros(2, N);  LW = L0;  s = zeros(1, N);
Wn = zeros(3, N);
for k = 1:N
    s(k) = ch3_phase(X(:, k), p0);
    [~, u0(:,k), a0] = ch3_col_dynamics(X(:, k), alpha, p0);
    [~, uW(:,k), aW] = ch3_col_dynamics(X(:, k), alpha, pW);
    L0(:,k) = a0.lam_drift + a0.lam_in * u0(:,k);
    LW(:,k) = aW.lam_drift + aW.lam_in * uW(:,k);
    [~, Wn(:,k)] = ch3_ub_wrench(X(1:7, k), s(k), pW);
end
names = {'stance hip', 'stance knee', 'swing hip', 'swing knee'};
fprintf('\nPART A  same gait, with vs without the upper-body wrench\n');
fprintf('  %-12s %10s %10s %12s\n', 'joint', 'peak|u0|', 'peak|uW|', 'RMS(du) [Nm]');
for j = 1:p0.nu
    fprintf('  %-12s %10.1f %10.1f %12.2f\n', names{j}, max(abs(u0(j,:))), ...
            max(abs(uW(j,:))), sqrt(mean((uW(j,:) - u0(j,:)).^2)));
end
fprintf('  GRF_z min: %.1f -> %.1f N,  max |Fx/Fz|: %.3f -> %.3f\n', ...
        min(L0(2,:)), min(LW(2,:)), max(abs(L0(1,:)./L0(2,:))), max(abs(LW(1,:)./LW(2,:))));

figure('Name', 'Upper body - part A');
subplot(3,1,1); plot(s, Wn.'); ylabel('wrench'); legend('F_x [N]','F_z [N]','M [Nm]');
title('Prescribed upper-body wrench on the torso'); grid on;
subplot(3,1,2); plot(s, u0.', '--', s, uW.', '-'); ylabel('u [Nm]');
title('Feedforward torque: dashed = none, solid = with wrench'); grid on;
subplot(3,1,3); plot(s, L0.', '--', s, LW.', '-'); ylabel('GRF [N]'); xlabel('phase s');
legend('F_x','F_z','F_x (W)','F_z (W)'); grid on;
drawnow;   % render now: MATLAB defers drawing while a long solve runs

%% ------------------------------------------------------------ PART B
fprintf('\nPART B  closed loop, %d steps, controller ''%s''\n', n_sim, ...
        char(getfield_or(p0, 'controller', 'default')));
x0   = X(:, 1);
out0 = ch3_simulate(x0, alpha, p0, n_sim);
outW = ch3_simulate(x0, alpha, pW, n_sim);
fprintf('  nominal : %d/%d steps ok  %s\n', out0.n_ok, n_sim, out0.reason);
fprintf('  wrench  : %d/%d steps ok  %s\n', outW.n_ok, n_sim, outW.reason);
report_speed('nominal', out0);  report_speed('wrench ', outW);

figure('Name', 'Upper body - part B');
lbl = {'x_{hip} [m]', 'z_{hip} [m]', 'q_t [rad]'};
for r = 1:3
    subplot(3,1,r);
    plot(out0.t, hip_row(out0, r), '--', outW.t, hip_row(outW, r), '-'); grid on;
    ylabel(lbl{r});
end
xlabel('t [s]'); legend('no wrench', 'with wrench');
drawnow;

%% ------------------------------------------------------------ PART C
if ~do_optim
    fprintf('\nPART C skipped: fmincon not found (Optimization Toolbox).\n');
    return;
end

% C1: known load, re-solve the gait (alpha, T, states) for it
fprintf('\nPART C1  re-optimizing the gait under the fixed wrench\n');
[z1, o1] = ch3_col_solve(pW, S.z);
fprintf('  exitflag %d, cost %.1f (nominal %.1f)\n', o1.exitflag, o1.fval, ch3_col_cost(S.z, p0));

% C2: crutch-like external support, solved with the gait
%     Fz >= 0 (a crutch can only push), |Fx| <= 60 N, no pure moment.
fprintf('\nPART C2  optimizing an external crutch-like wrench with the gait\n');
pC = ch3_ub_defaults(p0, 'bezier', 'optimize', true, 'deg', 4, ...
                     'lb', [-60; 0; 0], 'ub', [60; 200; 0], 'w', 500);
zC0 = ch3_col_pack(X, T, alpha, pC, [p0.theta_minus; p0.theta_plus]);
[zC, oC] = ch3_col_solve(pC, zC0);
pCe = oC.p;                                           % carries the solved beta
fprintf('  exitflag %d, total cost %.1f, of which wrench effort %.1f\n', ...
        oC.exitflag, oC.fval, ch3_ub_cost(pCe));
disp('  solved beta [Fx; Fz; M] (rows) :'); disp(pCe.ub.beta);
% Sweep pC.ub.w (e.g. 50, 500, 5000) to trace the torque-vs-support trade-off.

% ---------------------------------------------------------------- helpers
function v = getfield_or(p, f, d)
if isfield(p, f) && ~isempty(p.(f)), v = p.(f); else, v = d; end
end

function r = hip_row(out, k)
% x_hip = px, z_hip = -pz (pz is measured downward in q), q_t = q(3)
switch k
    case 1, r = out.x(1, :);
    case 2, r = -out.x(2, :);
    case 3, r = out.x(3, :);
end
end

function report_speed(tag, out)
if out.n_ok < 1, fprintf('  %s: no completed step\n', tag); return; end
dx = out.x(1, end) - out.x(1, 1);  dt = out.t(end) - out.t(1);
fprintf('  %s: mean speed %.3f m/s, q_t in [%.3f, %.3f] rad\n', tag, dx/dt, ...
        min(out.x(3, :)), max(out.x(3, :)));
end