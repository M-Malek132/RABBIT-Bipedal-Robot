function tf = ch4_is_stateful(p)
%CH4_IS_STATEFUL  Does p.controller carry internal state?
%
%   tf = ch4_is_stateful(p)
%
% The L1 laws are dynamical systems, not static feedback: they carry the state
% predictor, the two parameter estimates and the filter output described in
% ch4_l1_state. Everything else in Chapters 3 and 4 is a function of x alone.
%
% This is a one-line test, but it is asked in four places (ch4_step,
% ch4_simulate, ch4_forces, ch4_compare_controllers) and getting it wrong is
% silent -- a stateless path would simply never advance the adaptation and the
% L1 controller would degrade into the plain CLF-QP without complaining. So it
% lives in exactly one file.
%
% See also CH4_L1_STATE, CH4_CONTROL, CH4_STEP.

tf = any(strcmpi(p.controller, {'l1', 'l1_con'}));

end
