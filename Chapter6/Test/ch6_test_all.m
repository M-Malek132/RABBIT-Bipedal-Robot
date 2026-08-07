function pass = ch6_test_all()
%CH6_TEST_ALL  Run every Chapter-6 test.
%
%   pass = ch6_test_all()
%
% Order is the order to debug in -- each layer assumes the one below, so a
% failure high up is only worth reading once everything under it is green:
%
%   ch6_test_kin       the tracked points and rddot = Lf2r + LgLfr u, against
%                      the true flow. If this fails nothing else means anything:
%                      every barrier row would be enforcing a condition on
%                      something that is not an acceleration.
%   ch6_test_barrier   the GEOMETRY -- that inside O1 really does imply
%                      l_s <= l_max, that (6.19) really is tangent -- plus the
%                      derivative chain on a moving stone, and the identity
%                      between Section 6.1.1 and the rb = 2 ECBF of Section 5.2.
%   ch6_test_qp        the quadratic program: reduction to Chapter 3 at nb = 0,
%                      hard rows that are actually satisfied, the CLF as the
%                      only row that bends, and L_g g = 0 -- the measurement
%                      that motivates (6.1) in the first place.
%   ch6_test_foot3d    Section 6.3 on the swing-foot surrogate: the four
%                      barriers of Case 3 carried simultaneously.
%   ch6_test_sim       the stepping-stone harness: that it reproduces the
%                      reference gait when it should, that the barrier moves
%                      the placement when it should, and that the sampled hold
%                      converges.
%   ch6_test_library   the gait library map (6.22)-(6.24), on a synthetic
%                      library whose exact answer is known.
%
% Everything except ch6_test_sim runs in seconds. ch6_test_sim integrates
% several full walking steps with a QP at 1 kHz and takes a couple of minutes.
%
% See also CH6_MAIN, CH5_TEST_ALL, CH3_TEST_ALL.

fprintf('\n########## CHAPTER 6 TEST SUITE ##########\n');
t0 = tic;

suites = {@ch6_test_kin, @ch6_test_barrier, @ch6_test_qp, ...
          @ch6_test_foot3d, @ch6_test_sim, @ch6_test_library};
names  = {'kinematics', 'barriers', 'QP', 'foot3d', 'simulation', 'library'};

ok     = true;
failed = {};

for i = 1:numel(suites)
    try
        ok = suites{i}() && ok;
    catch ME
        ok = false;
        failed{end+1} = names{i};                               %#ok<AGROW>
        fprintf('\n  *** %s SUITE ERRORED: %s\n', upper(names{i}), ME.message);
        for k = 1:min(5, numel(ME.stack))
            fprintf('        at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
    end
end

fprintf('\n########## %s in %.1f s', upper(tf(ok)), toc(t0));
if isempty(failed)
    fprintf(' ##########\n\n');
else
    fprintf(' -- ERRORED: %s ##########\n\n', strjoin(failed, ', '));
end

if nargout > 0, pass = ok; end

end

function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
