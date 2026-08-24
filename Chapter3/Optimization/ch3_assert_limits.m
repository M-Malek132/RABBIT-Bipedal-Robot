function chk = ch3_assert_limits(z, p, what, logfile)
%CH3_ASSERT_LIMITS  Refuse to write a gait that violates its own enabled limits.
%
%   chk = ch3_assert_limits(z, p, what)
%   chk = ch3_assert_limits(z, p, what, logfile)
%
% Call this immediately BEFORE saving a deliverable gait.  It runs
% ch3_col_check_limits and throws if any inequality row that p itself declares
% enabled is violated beyond fmincon's constraint tolerance.
%
% WHY THIS IS AN ERROR RATHER THAN A WARNING.  A saved gait is not a diagnostic
% -- it is the input to everything downstream: ch3_report reads it, the Chapter
% 4/5/6 pipelines load it, and the marches in this folder warm-start from it.
% It carries the parameter struct that produced it, so a file whose p claims a
% constraint the z misses does not merely record a bad solve, it ASSERTS a
% false property, and every consumer inherits the assertion.  Warnings scroll
% past in a several-hour march; a refusal to write does not.
%
% WHAT IT IS NOT.  It is not a mesh check -- that is ch3_col_verify, and the
% two are independent.  Results/ch3_gait_forward_lean_tall.mat passed verify at
% 1.30e-05 while missing NEC3 by 0.92, which is precisely the gap this closes.
% Nor does it second-guess a solve that legitimately ran out of iterations: an
% infeasible fmincon result should stop the march at the stage that produced
% it, and this is the backstop for the case where it did not.
%
% Inputs
%   z, p    : the gait about to be written and the params to be written with it
%   what    : short description of the destination, used in the error message
%   logfile : optional; the failure report is appended here before throwing, so
%             an unattended march leaves the diagnosis on disk
%
% Output
%   chk     : the ch3_col_check_limits struct, for callers that want to record
%             it alongside the gait
%
% Throws ch3_assert_limits:violated.
%
% See also CH3_COL_CHECK_LIMITS, CH3_COL_VERIFY, CH3_UPGRADE_PARAMS.

if nargin < 3 || isempty(what), what = 'this gait'; end
if nargin < 4, logfile = ''; end

chk = ch3_col_check_limits(z, p);
if chk.ok, return; end

msg = sprintf(['Refusing to write %s.\n%s\n' ...
               'The parameter struct about to be saved beside this gait ' ...
               'declares those rows ENABLED, so writing it would assert a ' ...
               'property the trajectory does not have.\n' ...
               'Either re-solve with the gate on (ch3_impact_march for NEC3, ' ...
               'ch3_posture_march for the posture bands), or turn the gate ' ...
               'off in p so the file states what was actually enforced.'], ...
              what, chk.report);

if ~isempty(logfile)
    ch3_logln(logfile, msg);
end

error('ch3_assert_limits:violated', '%s', msg);

end
