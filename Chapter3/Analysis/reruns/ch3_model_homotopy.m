function ch3_model_homotopy(stage, seed_file, u_box)
%CH3_MODEL_HOMOTOPY  Carry a gait from the 30 kg robot to the 74 kg one.
%
%   ch3_model_homotopy()                   stage 'box', then stage 'free'
%   ch3_model_homotopy('box')              march lambda 0 -> 1 inside the box
%   ch3_model_homotopy('free')             continue from the last landed rung
%                                          with the torque gate off
%   ch3_model_homotopy(stage, seed_file, u_box)
%
% WHY. Every verified gait on today's 74 kg model descends from ONE warm start
% that jumped across the 2026-09-02 regeneration in a single step (Chapter 3,
% "the reference gait"): an old gait re-converged directly on the new
% dynamics, at 287 Nm, and then laddered down to 195 Nm and, from a 1.2 m/s
% gait, to a floor of 162.5 Nm. Nothing says where between the two robots a
% gait inside the declared 120 Nm stops existing, or why.
%
% A HOMOTOPY IN THE MODEL. The two robots have the same geometry, and M, V, G
% are linear in the link masses and inertias, so p.model_blend = lambda gives
% the robot with (1 - lambda) x the old and lambda x the new parameters,
% exactly (ch3_mvg). This walks lambda from 0 to 1, re-solving at every rung,
% warm-started from the rung below, WITH THE TORQUE BOX FIXED at u_box
% (default 120 Nm) and every other limit as the seed had it. The seed is
% Results/ch3_gait_full_constrained.mat: solved on the 30 kg model with all
% four Table 3.1 rows enforced at 120 Nm (peak 120.0 Nm, 0.975 m/s).
%
% Stage 'box' records, rung by rung, the gait that still fits (speed, stride,
% period, peak torque, impulse, contact) and stops where no gait does after
% the step has been halved three times: lambda*, the mass fraction at which
% the 120 Nm family ends, and the rows active just before it.
% Stage 'free' then continues from lambda* to 1 with the torque gate OFF, so
% the same family can be followed past the wall and the torque it NEEDS read
% off as the robot gets heavier -- the curve that says how far 120 Nm is from
% this family at 74 kg, instead of the single number "no gait below 195".
%
% Every landed rung is mesh-verified (ch3_col_verify) and saved with out.p,
% which records its model (model_sig, with the blend).
%
% Output: Results/reruns/model_homotopy/ (one .mat per rung, homotopy.log).
% Runtime: a warm-started N = 61 solve per rung, minutes to tens of minutes
% each; expect several hours for the whole walk. One MATLAB session at a time.
%
% See also CH3_MVG, CH3_COL_SOLVE, CH3_TORQUE_MARCH.

if nargin < 1 || isempty(stage), stage = {'box', 'free'}; end
if ischar(stage), stage = {stage}; end
ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
if nargin < 2 || isempty(seed_file)
    seed_file = fullfile(ROOT, 'Results', 'ch3_gait_full_constrained.mat');
end
if nargin < 3 || isempty(u_box), u_box = 120; end

SPD  = fullfile(ROOT, 'Results', 'reruns', 'model_homotopy');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'homotopy.log');
opts = struct('MaxFunctionEvaluations', 6e5, 'MaxIterations', 600);
DL0     = 0.05;                     % nominal lambda step
MIN_DL  = DL0 / 8;                  % three halvings, then stop

%% ------------------------------------------------------------ stage 'box'
if any(strcmp(stage, 'box'))
    S = load(seed_file);
    z = S.z(:);
    p = ch3_upgrade_params(S.p);
    p.scale_problem = true;
    p.model_blend   = 0;
    p.limits.u_max  = u_box;
    p.limits.enable.torque = true;
    logln(logf, '=== model homotopy, stage box | seed %s | box %.1f Nm | %s', ...
          seed_file, u_box, datestr(now)); %#ok<TNOW1,DATST>

    % Is the seed an orbit of the 30 kg robot (and, for contrast, of today's)?
    V0 = report(logf, p, z, 'seed @ lambda 0', NaN);
    p1 = p; p1.model_blend = 1;
    report(logf, p1, z, 'seed @ lambda 1', NaN);
    if ~V0.ok
        logln(logf, ['STOP: the seed is not a verified orbit of the 30 kg model, so ' ...
                     'there is nothing to continue. Pass another seed.']);
        return;
    end

    lam = 0; dl = DL0;
    while lam < 1 - 1e-12
        lt = min(1, lam + dl);
        f  = fullfile(SPD, sprintf('box_l%04.0f.mat', 1e4*lt));
        if exist(f, 'file')
            L = load(f);
            if L.landed
                z = L.z_try; p = L.p_try; lam = lt; dl = DL0;
                report(logf, p, z, sprintf('lambda %.4f (stored)', lt), NaN);
                continue;
            end
        else
            pk = p; pk.model_blend = lt;
            t0 = tic;
            [z_try, out] = ch3_col_solve(pk, z, opts);
            Vk = report(logf, out.p, z_try, sprintf('lambda %.4f', lt), toc(t0));
            landed = Vk.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
            p_try = out.p; %#ok<NASGU>
            save(f, 'z_try', 'p_try', 'out', 'Vk', 'landed', 'lt');
            if landed
                z = z_try; p = out.p; lam = lt; dl = DL0;
                continue;
            end
            logln(logf, '  rung %.4f missed (max|c| %.2e, max|ceq| %.2e, verify %d)', ...
                  lt, out.max_c, out.max_ceq, Vk.ok);
        end
        dl = dl / 2;
        if dl < MIN_DL
            logln(logf, ['=== WALL: no gait inside %.1f Nm beyond lambda = %.4f ' ...
                         '(%.1f kg); the last landed rung is %s'], u_box, lam, ...
                  30 + 44*lam, sprintf('box_l%04.0f.mat', 1e4*lam));
            break;
        end
    end
    if lam >= 1 - 1e-12
        logln(logf, '=== REACHED lambda = 1: a gait inside %.1f Nm exists on today''s robot', u_box);
    end
    lambda_star = lam; %#ok<NASGU>
    save(fullfile(SPD, 'box_result.mat'), 'lambda_star', 'z', 'p');
end

%% ----------------------------------------------------------- stage 'free'
if any(strcmp(stage, 'free'))
    R = fullfile(SPD, 'box_result.mat');
    if ~exist(R, 'file')
        logln(logf, 'stage free skipped: run stage box first'); return;
    end
    B = load(R);
    z = B.z; p = B.p; lam = B.lambda_star;
    if lam >= 1 - 1e-12
        logln(logf, 'stage free skipped: stage box reached lambda = 1'); return;
    end
    p.limits.enable.torque = false;             % follow the family past the wall
    logln(logf, '=== stage free | from lambda %.4f, torque gate off | %s', lam, ...
          datestr(now)); %#ok<TNOW1,DATST>
    dl = DL0;
    while lam < 1 - 1e-12
        lt = min(1, lam + dl);
        f  = fullfile(SPD, sprintf('free_l%04.0f.mat', 1e4*lt));
        if exist(f, 'file')
            L = load(f);
            if L.landed
                z = L.z_try; p = L.p_try; lam = lt; dl = DL0;
                report(logf, p, z, sprintf('free %.4f (stored)', lt), NaN);
                continue;
            end
        else
            pk = p; pk.model_blend = lt;
            t0 = tic;
            [z_try, out] = ch3_col_solve(pk, z, opts);
            Vk = report(logf, out.p, z_try, sprintf('free %.4f', lt), toc(t0));
            landed = Vk.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
            p_try = out.p;
            save(f, 'z_try', 'p_try', 'out', 'Vk', 'landed', 'lt');
            if landed
                z = z_try; p = p_try; lam = lt; dl = DL0;
                continue;
            end
        end
        dl = dl / 2;
        if dl < MIN_DL
            logln(logf, '=== stage free stops at lambda = %.4f', lam);
            break;
        end
    end
end
logln(logf, '=== MODEL_HOMOTOPY_DONE %s', datestr(now)); %#ok<TNOW1,DATST>
fprintf('MODEL_HOMOTOPY_DONE\n');
end

% ---------------------------------------------------------------------------
function V = report(logf, p, z, tag, secs)
E   = ch3_col_eval(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
[c, ceq] = ch3_col_constraints(z, p);
lam = [E.lam, E.lamm];
act = find(c(:).' > -1e-5);
blend = p.model_blend; if isempty(blend), blend = 1; end
logln(logf, ['%-22s %.1f kg | v %.4f | T %.4f | L %.4f | peak|u| %6.1f | impulse %5.2f | ' ...
             'mu %.3f | minFz %6.1f | max|ceq| %.1e | limits ok %d | verify %d (%.1e) | ' ...
             'active rows %s | %.0f s'], tag, 30 + 44*blend, E.L_step/E.T, E.T, E.L_step, ...
      max([abs(E.u(:)); abs(E.um(:))]), norm(E.impulse), max(abs(lam(1,:))./lam(2,:)), ...
      min(lam(2,:)), max(abs(ceq)), chk.ok, V.ok, V.max_dev, mat2str(act), secs);
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
