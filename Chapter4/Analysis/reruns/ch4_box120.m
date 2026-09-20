% Chapter 4 at the robot's DECLARED actuator limit: one 120 Nm box for every
% controller, every mass scale. This is the comparison sec:budget names as the
% chapter's remaining work.
%
% Resumable by design: one (controller, scale) run per iteration, each saved to
% its own file the moment it finishes, and an existing file is skipped. A
% libcurl abort therefore costs one run, not the sweep.

REPO = '/Users/mohammadmalek/Desktop/Research Files/RABBIT-Bipedal-Robot';
OUT  = fullfile(REPO, 'Results', 'box120');
LOG  = fullfile(REPO,'Results','box120','box120.log');
if ~exist(OUT, 'dir'), mkdir(OUT); end

U_BOX   = 120;
N_STEPS = 25;

jobs = { ...
    'robust', 'clfqp',      1.0; ...
    'robust', 'clfqp_con',  1.0; ...
    'robust', 'rclfqp_con', 1.0; ...
    'robust', 'clfqp',      1.5; ...
    'robust', 'clfqp_con',  1.5; ...
    'robust', 'rclfqp_con', 1.5; ...
    'robust', 'clfqp',      0.7; ...
    'robust', 'clfqp_con',  0.7; ...
    'robust', 'rclfqp_con', 0.7; ...
    'l1',     'l1',         1.0; ...
    'l1',     'l1_con',     1.0; ...
    'l1',     'l1',         0.7; ...
    'l1',     'l1_con',     0.7; ...
    'l1',     'l1',         1.5; ...
    'l1',     'l1_con',     1.5  ...
};

fid = fopen(LOG, 'a');
fprintf(fid, '=== box120 session start %s ===\n', datestr(now, 31));
fclose(fid);

[x0, alpha, p] = ch4_load_gait();

for j = 1:size(jobs, 1)
    preset = jobs{j,1};  name = jobs{j,2};  s = jobs{j,3};
    tag    = sprintf('%s_%s_s%03d', preset, name, round(s*100));
    f      = fullfile(OUT, [tag '.mat']);
    if exist(f, 'file')
        fid = fopen(LOG,'a'); fprintf(fid,'%-24s SKIP (done)\n', tag); fclose(fid);
        continue;
    end

    t0 = tic;
    opts = struct('controllers', {{name}}, 'scales', s, 'u_box', U_BOX, ...
                  'n_steps', N_STEPS, 'verbose', false, 'store_traj', false);
    try
        Cj = ch4_compare_controllers(x0, alpha, p, preset, opts);
        save(f, 'Cj', '-v7');
        fid = fopen(LOG,'a');
        fprintf(fid, '%-24s %6.1fs steps=%2d valid=%2d maxeta=%7.3f peak|u|=%7.1f minFz=%9.1f  %s\n', ...
                tag, toc(t0), Cj(1).steps_completed, Cj(1).valid_steps, ...
                Cj(1).max_eta, Cj(1).peak_torque, Cj(1).Fz_min, Cj(1).reason);
        fclose(fid);
    catch ME
        fid = fopen(LOG,'a');
        fprintf(fid, '%-24s %6.1fs ERROR %s | %s\n', tag, toc(t0), ME.identifier, ME.message);
        fclose(fid);
    end
end

fid = fopen(LOG,'a'); fprintf(fid, 'BOX120-ALL-DONE\n'); fclose(fid);
