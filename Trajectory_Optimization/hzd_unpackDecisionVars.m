function [CP, q0, dq0, T] = hzd_unpackDecisionVars(z, model, opt)
%HZD_UNPACKDECISIONVARS  Unpack  z = [CP_vec; q0; dq0; T]
%
%  CP_vec  : (n+1)*4 = opt.nCP_vars elements  (column-major flattening of CP)
%  q0      : nq = 7 elements
%  dq0     : nq = 7 elements
%  T       : 1 element
%
%  CP is recovered as  reshape(CP_vec, n+1, 4)
%    rows  = n+1 control points
%    cols  = 4   joints (q1 q2 q3 q4)

nCP = opt.nCP_vars;   % (n_bs+1)*ny
nq  = model.nq;

CP_vec = z(1 : nCP);
CP     = reshape(CP_vec, opt.n_bs+1, opt.ny);   % (n+1) x 4

idx = nCP;
q0  = z(idx+1 : idx+nq);   idx = idx + nq;
dq0 = z(idx+1 : idx+nq);   idx = idx + nq;
T   = z(idx+1);
end
