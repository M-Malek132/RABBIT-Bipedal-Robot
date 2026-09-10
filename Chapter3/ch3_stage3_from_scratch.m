%% ch3_stage3_from_scratch.m -- Stage 3 collocation solve, staged from a cold seed.
p = ch3_params();

U_TARGET = p.limits.u_max;        % captured ONCE, before anything mutates p
I_TARGET = p.limits.impulse_max;  % captured ONCE

z_seed = ch3_col_seed(p);
E_seed = ch3_col_eval(z_seed, p);
v_seed = E_seed.L_step / E_seed.T;
fprintf('seed natural speed = %.4f m/s (default v_des = %.2f)\n', v_seed, p.v_des);
p.v_des = v_seed * 1.1;
fprintf('using v_des = %.4f m/s\n', p.v_des);

staged = {'grf','friction','torque','impulse'};
fn = fieldnames(p.limits.enable);
for i = 1:numel(fn)
    p.limits.enable.(fn{i}) = ~ismember(fn{i}, staged);
end
for i = 1:numel(staged)
    p.limits.enable.(staged{i}) = false;
end

RESD = 'Results';

[z, out] = ch3_col_solve(p);
fprintf('stage bare: exitflag=%d ceq=%.3e c=%.3e\n', out.exitflag, out.max_ceq, out.max_c);
save(fullfile(RESD,'s3_00_bare.mat'), 'z','p','out');

p.limits.enable.grf = true;
[z, out] = ch3_col_solve(p, z);
fprintf('stage grf: exitflag=%d ceq=%.3e c=%.3e\n', out.exitflag, out.max_ceq, out.max_c);
save(fullfile(RESD,'s3_01_grf.mat'), 'z','p','out');

p.limits.enable.friction = true;
[z, out] = ch3_col_solve(p, z);
fprintf('stage friction: exitflag=%d ceq=%.3e c=%.3e\n', out.exitflag, out.max_ceq, out.max_c);
save(fullfile(RESD,'s3_02_friction.mat'), 'z','p','out');

%% torque ladder -- rungs computed from FIXED U_TARGET, never from p.limits.u_max
E = ch3_col_eval(z, p);
peak_u = max(max(abs(E.u(:))), max(abs(E.um(:))));
fprintf('natural peak |u| = %.1f Nm, target = %.1f Nm\n', peak_u, U_TARGET);
p.limits.enable.torque = true;
p.max_fun_evals = 450000;

if peak_u > U_TARGET
    start_u  = peak_u * 0.98;
    n_rungs  = max(3, ceil(log(start_u/U_TARGET) / log(1/0.75)) + 1);
    ratio_u  = (U_TARGET/start_u)^(1/(n_rungs-1));
    u_rungs  = start_u * ratio_u.^(0:n_rungs-1);
    u_rungs(end) = U_TARGET;   % force exact target on the last rung, avoid float drift
else
    u_rungs = U_TARGET;
end
fprintf('torque ladder (%d rungs): ', numel(u_rungs)); fprintf('%.1f ', u_rungs); fprintf('\n');

for k = 1:numel(u_rungs)
    p.limits.u_max = u_rungs(k);
    [z, out] = ch3_col_solve(p, z);
    fprintf('  torque[%d] u_max=%.2f: exitflag=%d ceq=%.3e c=%.3e\n', k, u_rungs(k), out.exitflag, out.max_ceq, out.max_c);
    save(fullfile(RESD, sprintf('s3_03_torque_%02d.mat', k)), 'z','p','out');
end

%% impulse ladder -- rungs computed from FIXED I_TARGET
E = ch3_col_eval(z, p);
peak_imp = norm(E.impulse);
fprintf('natural impulse = %.2f Ns, target = %.2f Ns\n', peak_imp, I_TARGET);
p.limits.enable.impulse = true;

if peak_imp > I_TARGET
    start_i  = peak_imp * 0.98;
    n_rungs  = max(3, ceil(log(start_i/I_TARGET) / log(1/0.75)) + 1);
    ratio_i  = (I_TARGET/start_i)^(1/(n_rungs-1));
    i_rungs  = start_i * ratio_i.^(0:n_rungs-1);
    i_rungs(end) = I_TARGET;
else
    i_rungs = I_TARGET;
end
fprintf('impulse ladder (%d rungs): ', numel(i_rungs)); fprintf('%.2f ', i_rungs); fprintf('\n');

for k = 1:numel(i_rungs)
    p.limits.impulse_max = i_rungs(k);
    [z, out] = ch3_col_solve(p, z);
    fprintf('  impulse[%d] max=%.2f: exitflag=%d ceq=%.3e c=%.3e\n', k, i_rungs(k), out.exitflag, out.max_ceq, out.max_c);
    save(fullfile(RESD, sprintf('s3_04_impulse_%02d.mat', k)), 'z','p','out');
end

verify = ch3_col_verify(z, p, false);
fprintf('FINAL: verify ok=%d max_dev=%.3e\n', verify.ok, verify.max_dev);
save(fullfile(RESD,'s3_final.mat'), 'z','p','out','verify');
