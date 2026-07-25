function [z2, p2] = ch3_col_remesh(z, p, N_new)
%CH3_COL_REMESH  Move a collocation solution onto a finer (or coarser) mesh.
%
%   [z2, p2] = ch3_col_remesh(z, p, N_new)
%
% Mesh refinement is the standard cure for a spurious discrete solution: solve
% cheaply on a coarse mesh, then re-solve on a finer one warm-started from the
% result, until ch3_col_verify stops complaining.  Going straight to a fine
% mesh from a cold start is both slower per iteration and more likely to stall,
% because the coarse solve does the cheap structural work first.
%
% The node states are re-interpolated with pchip in the normalized step time
% tau = t/T, so the same routine works whether or not T changed. alpha and T
% carry over untouched -- alpha is mesh-independent by construction, which is
% precisely why parametrizing the gait by Bezier coefficients rather than by
% the node states is worth doing.
%
% Inputs
%   z, p  : solution and parameters on the current mesh
%   N_new : node count for the new mesh
%
% Outputs
%   z2 : packed decision vector on the new mesh
%   p2 : p with N_nodes updated (pass this to ch3_col_solve)
%
% See also CH3_COL_VERIFY, CH3_COL_SOLVE.

[X, T, alpha] = ch3_col_unpack(z, p);
N = size(X, 2);

tau_old = linspace(0, 1, N);
tau_new = linspace(0, 1, N_new);

X2 = interp1(tau_old(:), X.', tau_new(:), 'pchip').';

p2 = p;
p2.N_nodes = N_new;

z2 = ch3_col_pack(X2, T, alpha, p2);

end
