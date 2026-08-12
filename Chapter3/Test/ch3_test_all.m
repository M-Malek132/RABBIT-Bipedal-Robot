function ch3_test_all()
%CH3_TEST_ALL  Run every Chapter-3 verification suite.
%
%   ch3_test_all
%
% Each suite checks its stage against an INDEPENDENT reference rather than
% against itself:
%
%   ch3_test_params       ch3_upgrade_params, which every analysis and every
%                         warm start routes a saved p through: missing fields
%                         refilled from the defaults, present fields preserved,
%                         and the one field deliberately NOT preserved.
%   ch3_test_model        f + g*u vs the KKT constrained dynamics; impact
%                         against a hand-built momentum balance; guard
%                         gradient against finite differences.
%   ch3_test_vc           Bezier derivatives against finite differences; the
%                         Bezier/B-spline equivalence; and ydd = Lf^2y + LgLfy u
%                         against the derivative of the TRUE flow.
%   ch3_test_control      the (C)ARE residuals; Vdot affine in mu; the stage-7
%                         closed form against quadprog; and each controller
%                         against its OWN Lyapunov certificate.
%   ch3_test_collocation  pack/unpack; seed quality; degrees of freedom; and
%                         the Hermite-Simpson order of accuracy under mesh
%                         refinement.
%   ch3_test_simulation   the zero-order-hold integrator against the continuous
%                         rollout of the same controller, checking that the
%                         sampling error is first order in the period.
%   ch3_test_hzd          the Section 6.3.4 constraint set: the reduction to
%                         the scalar zero dynamics against the full 14-state
%                         model, the conserved (1/2)sigma^2 + V_zero along a
%                         real rollout, and delta_zero^2 against the Poincare
%                         spectral radius from 26 independent simulations.
%
% Runs in a few minutes. Anything that fails prints the measured error and its
% tolerance, so the size of the discrepancy is visible rather than just a name.

suites = {@ch3_test_params, @ch3_test_model, @ch3_test_vc, @ch3_test_control, ...
          @ch3_test_collocation, @ch3_test_simulation, @ch3_test_hzd};
names  = {'params', 'model', 'virtual constraints', 'controllers', ...
          'collocation', 'simulation', 'hybrid zero dynamics'};

fprintf('\n############ CHAPTER 3 TEST SUITE ############\n');
t0 = tic;
failed = {};

for i = 1:numel(suites)
    try
        suites{i}();
    catch ME
        failed{end+1} = names{i}; %#ok<AGROW>
        fprintf('\n  *** %s SUITE ERRORED: %s\n', upper(names{i}), ME.message);
        for k = 1:numel(ME.stack)
            fprintf('        at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
    end
end

fprintf('\n############ done in %.1f s', toc(t0));
if isempty(failed)
    fprintf(' ############\n\n');
else
    fprintf(' -- ERRORED: %s ############\n\n', strjoin(failed, ', '));
end

end
