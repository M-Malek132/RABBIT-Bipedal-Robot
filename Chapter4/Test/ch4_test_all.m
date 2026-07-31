function ch4_test_all()
%CH4_TEST_ALL  Run every Chapter-4 test.
%
%   ch4_test_all()
%
% Order is deliberate and is the order to debug in: the model layer first
% (nothing above it can be right if the true/nominal split is wrong), then the
% robust controller, then the adaptive one.
%
%   ch4_test_model   true vs nominal dynamics, and the Delta terms of (4.4)
%   ch4_test_rclf    the robust CLF-QP of Section 4.1 and its guarantee
%   ch4_test_l1      the L1 adaptive controller of Section 4.2
%
% See also CH3_TEST_ALL, CH4_MAIN.

fprintf('\n########## CHAPTER 4 TEST SUITE ##########\n');
t0 = tic;

ch4_test_model();
ch4_test_rclf();
ch4_test_l1();

fprintf('########## done in %.1f s ##########\n\n', toc(t0));

end
