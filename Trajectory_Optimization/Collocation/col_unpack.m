function [X, U, Lam, coeffs, T] = col_unpack(z, p)
%COL_UNPACK  Split the collocation decision vector into its blocks.
%
%   [X, U, Lam, coeffs, T] = col_unpack(z, p)
%
% Layout of z (see col_pack):
%   X      2*nq x N   node states [q; dq]
%   U      nu   x N   node joint torques
%   Lam    2    x N   node stance-foot contact forces [Fx; Fz]
%   coeffs nu x n_coeffs   B-spline virtual-constraint parameters
%   T      scalar     step duration
%
% p must carry p.N (node count), p.nq, p.nu, p.n_coeffs.

    N  = p.N;  nq = p.nq;  nu = p.nu;  nc = p.n_coeffs;
    nx = 2*nq;

    iX = nx*N;   iU = nu*N;   iL = 2*N;   iC = nu*nc;

    X      = reshape(z(1:iX),                 nx, N);
    U      = reshape(z(iX + (1:iU)),          nu, N);
    Lam    = reshape(z(iX + iU + (1:iL)),      2, N);
    coeffs = reshape(z(iX + iU + iL + (1:iC)), nu, nc);
    T      = z(iX + iU + iL + iC + 1);
end
