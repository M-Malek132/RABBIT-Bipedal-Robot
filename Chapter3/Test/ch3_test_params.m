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
% ...with exactly two documented exceptions, checked below.
%
%   1-3.  missing fields, including nested limits, come from the defaults.
%   3b.   EXCEPT p.limits.enable, where a missing gate is filled in as FALSE.
%   4-5.  present fields win, and are not overwritten by the defaults.
%   6.    p.checkpoint_file is cleared when this machine cannot honour it.
%
% WHY 3b EXISTS.  Filling a missing gate from ch3_params makes the loader's
% answer depend on a default that can change later -- and it did.  b64160e
% added the six NIC/NEC gates defaulting to false; e7e101b flipped every
% default to true; from that commit on, eight gaits solved before the gates
% existed loaded claiming to enforce them.  Results/ch3_gait_forward_lean_tall
% .mat is the sharp case: it verifies as a real trajectory at 1.30e-05 and
% misses NEC3 by 0.92, and it is the documented warm-start seed for
% ch3_lean_tall_march.  Nothing was written and nothing re-solved to cause
% that; only a default in another file moved.
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

%% 3b. ...BUT REFILLED AS FALSE, which is the opposite of every other field.
% A gate absent from a saved struct is evidence the constraint was not enforced
% when that gait was solved, not an opinion the result forgot to record. Taking
% it from the defaults instead is what silently switched six constraints on
% underneath eight stored gaits when e7e101b flipped every default to true --
% Results/ch3_gait_forward_lean_tall.mat then loaded claiming to enforce NEC3
% while missing it by 0.92. These cases pin the rule so a future default flip
% cannot reach back into results already on disk.
p = d;  p.limits.enable = rmfield(p.limits.enable, en{1});
q = ch3_upgrade_params(p);
pass = ok('missing enable flag refilled as FALSE', ...
          q.limits.enable.(en{1}) == false, pass);

% The historical shape exactly: the six-field enable struct that shipped before
% b64160e added the NIC/NEC gates.
p = d;
p.limits.enable = struct('torque', false, 'impulse', false, 'friction', false, ...
                         'grf', false, 'clearance', true, 'height', true);
q = ch3_upgrade_params(p);
pass = ok('pre-NEC struct gets no phantom gates', ...
          ~q.limits.enable.impact && ~q.limits.enable.hzd && ...
          ~q.limits.enable.swing_clear && ~q.limits.enable.liftoff && ...
          ~q.limits.enable.phase_mono && ~q.limits.enable.decoupling, pass);
pass = ok('...while its own six flags survive', ...
          q.limits.enable.clearance && q.limits.enable.height && ...
          ~q.limits.enable.torque && ~q.limits.enable.grf, pass);

% An explicit true is still honoured -- the rule is "absent means off", not
% "everything off".
p = d;
p.limits.enable = struct('impact', true);
q = ch3_upgrade_params(p);
pass = ok('explicit true still wins', ...
          q.limits.enable.impact && ~q.limits.enable.torque, pass);

% A struct predating gating altogether.
p = d;  p.limits = rmfield(p.limits, 'enable');
q = ch3_upgrade_params(p);
pass = ok('enable absent entirely -> all gates off', ...
          isfield(q.limits, 'enable') && ...
          ~any(struct2array_local(q.limits.enable)), pass);

% And a FRESH ch3_params must round-trip untouched: the rule applies to what a
% saved struct omits, and a fresh one omits nothing.
q = ch3_upgrade_params(d);
pass = ok('fresh ch3_params round-trips unchanged', ...
          isequal(q.limits.enable, d.limits.enable), pass);

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

function v = struct2array_local(s)
% struct2array ships with a toolbox this suite does not require.
c = struct2cell(s);
v = [c{:}];
end
