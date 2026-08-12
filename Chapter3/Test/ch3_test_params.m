function ch3_test_params()
%CH3_TEST_PARAMS  Stage-0 verification: loading a saved parameter struct.
%
% Results are saved with the p that produced them, so every analysis and every
% warm start starts by feeding an OLD struct to ch3_upgrade_params.  That makes
% this function load-bearing for the whole chapter, and it has two jobs pulling
% in opposite directions:
%
%   1. FILL IN what a stale struct is missing, so analysis code can read fields
%      added after the result was written without erroring.
%   2. PRESERVE what the stale struct has, so a result is still analysed under
%      the settings that produced it.
%
% ...with exactly one documented exception, checked below.
%
%   1-3. missing fields, including nested limits, come from the defaults.
%   4-5. present fields win, and are not overwritten by the defaults.
%   6.   p.checkpoint_file is cleared when this machine cannot honour it.
%
% WHY 6 EXISTS.  A .mat saved on Windows carries an absolute
% 'C:\Users\...\Results\ch3_*_ckpt.mat'.  On macOS and Linux a backslash is a
% LEGAL FILENAME CHARACTER, so ch3_col_solve's checkpoint save does not fail --
% it succeeds, creating a single file whose entire Windows path is its name, in
% whatever the working directory happened to be.  Two of those were committed to
% this repo before anyone noticed, from two different directories.  The
% try/catch around that save only catches an UNWRITABLE destination; this one is
% writable and simply wrong, so ch3_upgrade_params has to catch it instead.
%
% See also CH3_UPGRADE_PARAMS, CH3_PARAMS, CH3_COL_SOLVE.

fprintf('\n=== ch3_test_params ===\n');
pass = true;
d    = ch3_params();

%% 1-3. missing fields are filled in from the defaults
p = rmfield(d, 'v_des');
q = ch3_upgrade_params(p);
pass = ok('missing scalar refilled', isfield(q, 'v_des') && q.v_des == d.v_des, pass);

p = d;  p.limits = rmfield(p.limits, 'mu_s');
q = ch3_upgrade_params(p);
pass = ok('missing nested limit refilled', ...
          isfield(q.limits, 'mu_s') && q.limits.mu_s == d.limits.mu_s, pass);

en = fieldnames(d.limits.enable);
p = d;  p.limits.enable = rmfield(p.limits.enable, en{1});
q = ch3_upgrade_params(p);
pass = ok('missing enable flag refilled', ...
          isequal(sort(fieldnames(q.limits.enable)), ...
                  sort(fieldnames(d.limits.enable))), pass);

%% 4-5. present fields are preserved, NOT overwritten by the defaults
p = d;  p.v_des = 0.4242;  p.N_nodes = 17;
q = ch3_upgrade_params(p);
pass = ok('saved scalars still win', ...
          q.v_des == 0.4242 && q.N_nodes == 17, pass);

p = d;  p.limits.u_max = 999;
q = ch3_upgrade_params(p);
pass = ok('saved nested limit still wins', q.limits.u_max == 999, pass);

%% 6. checkpoint_file is the one field deliberately not preserved
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
good = fullfile(root, 'Results', 'ch3_unit_ckpt.mat');   % dir exists; file need not
gone = fullfile(root, 'Results', 'no_such_dir_xyz', 'a.mat');

% Every expectation below holds on BOTH platforms, so the suite is not
% conditioned on ispc.  The two backslash cases name a directory that does not
% exist anywhere: on POSIX they are cleared by the backslash rule, on Windows by
% the missing-directory rule.  Either way the answer is ''.
cases = { ...
    'windows absolute cleared', 'C:\no_such_dir_xyz\Results\ch3_posture_ckpt.mat', ''; ...
    'windows relative cleared', 'no_such_dir_xyz\ch3_x_ckpt.mat',                   ''; ...
    'stale directory cleared',  gone,                                               ''; ...
    'empty stays empty',        '',                                                 ''; ...
    'usable path preserved',    good,                                             good};

for k = 1:size(cases, 1)
    p = d;  p.checkpoint_file = cases{k,2};
    q = ch3_upgrade_params(p);
    pass = ok(cases{k,1}, strcmp(q.checkpoint_file, cases{k,3}), pass);
end

fprintf('--- ch3_test_params: %s ---\n\n', tf(pass));
end

% ---------------------------------------------------------------------------
function pass = ok(name, cond, pass_in)
pass = pass_in && cond;
fprintf('  [%s] %-34s\n', tf(cond), name);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
