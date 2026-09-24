function chk = ch3_col_check_limits(z, p, tol)
%CH3_COL_CHECK_LIMITS  Does this gait satisfy the limits its own p enables?
%
%   chk = ch3_col_check_limits(z, p)
%   chk = ch3_col_check_limits(z, p, tol)
%
% THE QUESTION CH3_COL_VERIFY DOES NOT ASK.
%
% ch3_col_verify asks whether the node states lie on a real trajectory of the
% closed-loop ODE.  That is a question about the MESH, and a gait can pass it
% perfectly while sitting far outside a constraint its own parameter struct
% declares enabled -- the two are independent failures.  Measured on this
% project: Results/ch3_gait_forward_lean_tall.mat verified at 1.30e-05 (ok = 1)
% with max|ceq| = 1.04e-07, and simultaneously missed NEC3, the impact friction
% cone, by 0.92.  Every number on that file looked healthy because nothing was
% asking the second question.
%
% This asks it.  Evaluate the inequality block at z under p, and report every
% row that is positive beyond the tolerance fmincon was told to hold.  A gait
% that fails this is not "slightly off" -- it is a minimizer of a DIFFERENT
% problem than the one p describes, and warm-starting from it begins infeasible
% in whatever direction the missing row constrains.
%
% WHY THE GATES MATTER AND DISABLED ROWS DO NOT.  ch3_col_constraints holds a
% disabled row at -1 rather than removing it, so c always has length 19 and the
% row indices below are stable regardless of what is on.  A disabled row is
% therefore trivially satisfied and cannot appear here: this function only ever
% complains about a limit the params themselves claim to be enforcing, which is
% exactly the self-consistency property a stored gait should have.
%
% Inputs
%   z   : decision vector
%   p   : parameter struct (upgraded here, so a loaded p is fine)
%   tol : violation threshold.  Default p.con_tol if present, else 1e-6 --
%         fmincon's ConstraintTolerance in ch3_col_solve, which is the level a
%         converged solve is entitled to be believed at.
%
% Output
%   chk : struct
%       .ok        true if no ENABLED inequality is violated beyond tol
%       .max_c     max(c) over all rows (disabled rows sit at -1)
%       .rows      indices of the violated rows
%       .names     cellstr describing them, one per .rows entry
%       .values    c(.rows)
%       .tol       the tolerance used
%       .max_ceq   max|ceq|, reported for context (equalities are not gated)
%       .report    printable multi-line summary, '' when ok
%
% See also CH3_COL_CONSTRAINTS, CH3_COL_VERIFY, CH3_UPGRADE_PARAMS.

p = ch3_upgrade_params(p);

if nargin < 3 || isempty(tol)
    if isfield(p, 'con_tol') && ~isempty(p.con_tol)
        tol = p.con_tol;
    else
        tol = 1e-6;
    end
end

% Row labels track the table in ch3_col_constraints' header. Kept here rather
% than returned by that function so the constraint evaluation stays a hot path
% with no cellstr construction in it -- fmincon calls it once per finite
% difference, i.e. thousands of times per iteration.
names = { ...
    'mid-step swing-foot clearance          [clearance]', ...
    'no ground penetration at interior nodes[always]', ...
    'step-length floor                      [always]', ...
    'peak torque  |u| <= u_max              [torque]', ...
    'friction cone |Fx| <= mu_s Fz     NIC2 [friction]', ...
    'minimum normal force Fz >= Fz_min NIC1 [grf]', ...
    'impact impulse ||I|| <= impulse_max    [impulse]', ...
    'hip-height band                        [height]', ...
    'swing foot strictly above ground  NIC3 [swing_clear]', ...
    'transversal strike (descending at N)NIC3[swing_clear]', ...
    'post-impact swing-leg lift-off    NEC2 [liftoff]', ...
    'impact impulse compressive Iz>=0  NEC3 [impact]', ...
    'impact impulse inside friction cone NEC3[impact]', ...
    'existence of the fixed point      NEC4 [hzd]', ...
    'stability of the fixed point      NEC5 [hzd]', ...
    'theta strictly monotonic          HH6  [phase_mono]', ...
    'decoupling matrix invertible on Z HH2  [decoupling]', ...
    'swing-foot height ceiling              [clearance_max]', ...
    'step-length ceiling                    [step_len_max]'};

[c, ceq] = ch3_col_constraints(z, p);

rows = find(c(:).' > tol);

% Built field by field, NOT by struct(...) with a cell argument.  struct() with
% a 1x0 cell for .names would return a 1x0 STRUCT ARRAY rather than a scalar
% struct with an empty cell field -- so the ok case, the common one, would come
% back as something with no fields to read.
chk         = struct();
chk.ok      = isempty(rows);
chk.max_c   = max(c);
chk.rows    = rows;
chk.names   = names(rows);
chk.values  = c(rows);
chk.tol     = tol;
chk.max_ceq = max(abs(ceq));
chk.report  = '';

if chk.ok, return; end

lines = cell(1, numel(rows) + 1);
lines{1} = sprintf(['constraint-consistency FAILED: %d enabled ' ...
                    'inequality row(s) violated beyond %.1e'], numel(rows), tol);
for i = 1:numel(rows)
    lines{i+1} = sprintf('    c(%2d) = %+.4e   %s', ...
                         rows(i), c(rows(i)), names{rows(i)});
end
chk.report = strjoin(lines, newline);

end
