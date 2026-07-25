function z = ch3_col_pack(X, T, alpha, p)
%CH3_COL_PACK  Assemble the collocation decision vector.
%
%   z = ch3_col_pack(X, T, alpha, p)
%
% Layout:
%
%       z = [ X(:) ; T ; alpha(:) ]
%           \_ 14N _/   1   \_ ny*n_ctrl _/
%
%   X     : 14 x N  state at every collocation node
%   T     : scalar  step duration (free -- the gait finds its own timing)
%   alpha : ny x n_ctrl  the Bezier coefficients, which are what stage 3 is
%           actually solving for; the node states are there to make the
%           dynamics algebraic rather than an ODE inside the optimizer.
%
% See also CH3_COL_UNPACK.

z = [X(:); T; alpha(:)];

if numel(z) ~= p.nx*size(X,2) + 1 + p.ny*p.n_ctrl
    error('ch3_col_pack:size', 'Packed vector has unexpected length.');
end

end
