function [u, info] = ch5_control(x, p, e)
%CH5_CONTROL  Evaluate the selected Chapter-5 controller at a state.
%
%   [u, info] = ch5_control(x, p)
%   [u, info] = ch5_control(x, p, e)     with a prebuilt ch5_ecbf_gain
%
% The single entry point every simulation and test uses. All four controllers
% share the same two upstream objects -- the input-output linearization of the
% tracking outputs, and the Lie stack of the constraint -- and differ only in
% which rows they put in the quadratic program:
%
%   'clfqp'          CLF row only                       Fig 5.3a / 5.4a
%   'cbfclfqp'       CLF + reciprocal barrier (5.7)     Section 5.1
%   'cbfclfqp_viol'  the same via VIOL (5.15)           Remark 5.4
%   'ecbfclfqp'      CLF + exponential barrier (5.31)   Section 5.2
%   'ecbfclfqp_viol' the same with mu_b explicit        Remark 5.4 again
%
% THE BARRIER IS EVALUATED EVEN FOR 'clfqp'. It costs almost nothing and it
% means every run in this chapter carries the same diagnostic record --
% h(t), eta_b(t), and whether the row would have been active -- whether or not
% the controller was listening to it. That is what makes the baseline's
% violation measurable rather than merely visible on a plot.
%
% Inputs
%   x : nx x 1 state
%   p : parameter struct
%   e : optional ECBF gain struct; rebuilt from p if omitted
%
% Outputs
%   u    : nu x 1 control
%   info : struct .io .b .qp .mu .u .h .controller
%
% See also CH5_IO_LIN, CH5_BARRIER, CH5_CTRL_ECBF_CLF_QP, CH5_ODE_RHS.

x  = x(:);
io = ch5_io_lin(x, p);
b  = ch5_barrier(x, p);

switch lower(p.controller)

    case 'clfqp'
        [mu, u, qp] = ch5_ctrl_clf_qp(io, p);

    case 'cbfclfqp'
        [mu, u, qp] = ch5_ctrl_cbf_clf_qp(io, b, p, false);

    case 'cbfclfqp_viol'
        [mu, u, qp] = ch5_ctrl_cbf_clf_qp(io, b, p, true);

    case {'ecbfclfqp', 'ecbfclfqp_viol'}
        if nargin < 3 || isempty(e)
            e = ch5_ecbf_gain(p, b.rb);
        end
        [mu, u, qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, p, ...
                                           strcmpi(p.controller, 'ecbfclfqp_viol'));

    otherwise
        error('ch5_control:unknownController', ...
              ['Unknown p.controller "%s" (expected clfqp | cbfclfqp | ' ...
               'cbfclfqp_viol | ecbfclfqp | ecbfclfqp_viol).'], p.controller);
end

info = struct('io', io, 'b', b, 'qp', qp, 'mu', mu, 'u', u, 'h', b.h, ...
              'controller', lower(p.controller));

end
