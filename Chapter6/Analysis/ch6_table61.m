function T = ch6_table61(p, alpha_nom, lib)
%CH6_TABLE61  The three-controller comparison of (6.27) and Table 6.1.
%
%   T = ch6_table61(p, alpha_nom)        controllers I and II only
%   T = ch6_table61(p, alpha_nom, lib)   all three
%
% For each step-length range, generate p.mc.n_trials random terrains of
% p.mc.n_stones stones each and run all three controllers of (6.27) on the SAME
% terrains:
%
%   I   : Gait Library          -- interpolated gait, NO barrier rows
%   II  : CBF                   -- one nominal gait, barrier rows
%   III : CBF and Gait Library  -- both
%
% and report the percentage of trials in which the robot crossed the whole
% terrain without violating foot placement, ground contact, friction or input
% constraints. That is the thesis's own success criterion, applied by
% ch6_simulate.
%
% ================================================== THE SAME TERRAINS, ON PURPOSE
% All three controllers see identical terrain for a given (range, trial). The
% table's content is a DIFFERENCE between controllers, and independent random
% draws would add sampling noise to exactly the quantity being measured -- at
% 20 trials that noise is comparable to the differences in the middle rows.
% Seeding from p.mc.seed and the trial index gives paired samples for free and
% makes any individual failing trial re-runnable on its own.
%
% ==================================================== WHAT IS AND IS NOT RUN HERE
% The thesis runs 100 problem sets x 10 stones for 7 ranges and 3 controllers:
% 21000 simulated walking steps, each an event-terminated integration with a QP
% at 1 kHz. That is days of compute. p.mc.n_trials and p.mc.ranges are therefore
% knobs, and this function PRINTS the counts it used at the top of the table so
% a reduced run can never be mistaken for the full one.
%
% A reduced trial count widens the confidence interval; it does not bias the
% estimate. The binomial standard error at n trials is at most 0.5/sqrt(n) --
% 11 percentage points at n = 20, 5 at n = 100 -- so it is printed alongside
% each cell. Differences smaller than that are not differences.
%
% Inputs
%   p         : parameter struct (uses p.mc)
%   alpha_nom : the one nominal gait, for controllers I and II
%   lib       : gait library struct; omit to skip controller III
%
% Output
%   T : struct
%       .ranges   1 x nr cell
%       .pct      nr x 3 success percentage
%       .se       nr x 3 binomial standard error, percentage points
%       .n_trials .n_stones .stone_sz
%       .detail   nr x 3 cell of per-trial reasons
%
% See also CH6_SIMULATE, CH6_TERRAIN, CH6_LIB_ALPHA, CH6_REPORT.

if nargin < 3, lib = []; end

nr = numel(p.mc.ranges);
have_lib = ~isempty(lib);
nc = 2 + have_lib;

names = {'Gait Library', 'CBF', 'CBF & Gait Library'};
names = names(1:nc);
if ~have_lib, names = {'Gait Library', 'CBF'}; end

pct    = zeros(nr, nc);
se     = zeros(nr, nc);
detail = cell(nr, nc);

fprintf('\n=============== TABLE 6.1 (as run here) ===============\n');
fprintf(' %d trials x %d stones, stone size %.0f cm, seed %d\n', ...
        p.mc.n_trials, p.mc.n_stones, 100*p.mc.stone_sz, p.mc.seed);
fprintf(' thesis runs 100 x 10 over 7 ranges; this is a reduced grid\n');
fprintf(' controller: %s barrier, poles (%.4g, %.4g), torque box %s\n', ...
        p.cbf.form, p.cbf.gamma_b, p.cbf.gamma, onoff(p.limits.enable.torque));

for r = 1:nr
    rng_r = p.mc.ranges{r};

    for c = 1:nc
        [q, gait] = configure(p, c, alpha_nom, lib, have_lib);

        n_ok = 0;
        why  = cell(1, p.mc.n_trials);
        for tr = 1:p.mc.n_trials
            seed = p.mc.seed + 1000*r + tr;
            terr = ch6_terrain(p.mc.n_stones, rng_r, p.mc.stone_sz, seed);

            s = ch6_simulate(nominal_x0(), gait, q, terr);

            n_ok  = n_ok + s.success;
            why{tr} = s.reason;
        end

        pct(r,c)    = 100 * n_ok / p.mc.n_trials;
        se(r,c)     = 100 * sqrt(max(pct(r,c)/100 * (1-pct(r,c)/100), 0) ...
                                 / p.mc.n_trials);
        detail{r,c} = why;
        fprintf('   [%.2f:%.2f] %-20s %5.0f%% +- %.0f\n', ...
                rng_r(1), rng_r(2), names{c}, pct(r,c), se(r,c));
    end
end

%% ------------------------------------------------------------------- table
fprintf('\n %-18s', 'Step Length (cm)');
for c = 1:nc, fprintf(' %-20s', names{c}); end
fprintf('\n');
for r = 1:nr
    fprintf(' [%2.0f:%2.0f]%10s', 100*p.mc.ranges{r}(1), 100*p.mc.ranges{r}(2), '');
    for c = 1:nc
        fprintf(' %3.0f%% (+-%2.0f)%9s', pct(r,c), se(r,c), '');
    end
    fprintf('\n');
end
fprintf('=======================================================\n');

T = struct('ranges', {p.mc.ranges}, 'names', {names}, ...
           'pct', pct, 'se', se, 'detail', {detail}, ...
           'n_trials', p.mc.n_trials, 'n_stones', p.mc.n_stones, ...
           'stone_sz', p.mc.stone_sz, 'seed', p.mc.seed);

end

% ---------------------------------------------------------------------------
function [q, gait] = configure(p, c, alpha_nom, lib, have_lib)
%CONFIGURE  The three controllers of (6.27), differing only where they should.
%
% Same integrator, same sampling rate, same torque box, same limits. The ONLY
% differences are the barrier rows and where alpha comes from -- otherwise the
% columns would not be comparable and the table would be measuring the harness.
q = p;
switch c
    case 1                      % I: gait library, no CBF
        q.controller  = 'clfqp_con';
        q.cbf.problem = 'none';
        gait = pick_lib(lib, alpha_nom, have_lib);
    case 2                      % II: CBF, one nominal gait
        q.controller  = 'cbf_clf_qp';
        q.cbf.problem = 'stones';
        gait = alpha_nom;
    case 3                      % III: both
        q.controller  = 'cbf_clf_qp';
        q.cbf.problem = 'stones';
        gait = lib;
end
end

function g = pick_lib(lib, alpha_nom, have_lib)
if have_lib, g = lib; else, g = alpha_nom; end
end

% ---------------------------------------------------------------------------
function x0 = nominal_x0()
%NOMINAL_X0  Every trial starts from the reference gait's fixed point.
%
% The same start state for every controller and every trial, so a difference
% between columns is a difference between controllers. Cached because a table
% runs it hundreds of times and it is a file read.
persistent X0
if isempty(X0)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    S = load(fullfile(root, 'Results', 'ch3_reference_gait.mat'));
    X  = ch3_col_unpack(S.z_opt, S.p);
    X0 = X(:,1);
end
x0 = X0;
end

function s = onoff(b)
if b, s = 'on'; else, s = 'off'; end
end
