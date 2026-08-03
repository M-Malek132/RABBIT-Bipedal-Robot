function pass = ch5_test_all()
%CH5_TEST_ALL  Run every Chapter-5 test.
%
%   pass = ch5_test_all()
%
% Order is deliberate and is the order to debug in. Each layer assumes the one
% below it, so a failure high up is only worth reading once everything under it
% is green:
%
%   ch5_test_model    the two plants, and the RELATIVE DEGREES the chapter
%                     rests on -- measured from the model, not declared
%   ch5_test_barrier  the Lie stack eta_b of (5.10), including the committed
%                     symbolic code, checked against finite differences taken
%                     along an actual drift trajectory
%   ch5_test_ecbf     Section 5.2's theory: pole placement, Remark 5.6,
%                     Proposition 5.1, Corollary 5.2, and the guarantee itself
%   ch5_test_qp       the quadratic programs: Remark 5.4's equivalences, and
%                     the high-relative-degree failure that motivates the
%                     whole chapter
%
% See also CH5_MAIN, CH4_TEST_ALL, CH3_TEST_ALL.

fprintf('\n########## CHAPTER 5 TEST SUITE ##########\n');
t0 = tic;

ok = true;
ok = ch5_test_model()   && ok;
ok = ch5_test_barrier() && ok;
ok = ch5_test_ecbf()    && ok;
ok = ch5_test_qp()      && ok;

fprintf('########## %s in %.1f s ##########\n\n', ...
        upper(tf(ok)), toc(t0));

if nargout > 0, pass = ok; end

end

% ---------------------------------------------------------------------------
function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
