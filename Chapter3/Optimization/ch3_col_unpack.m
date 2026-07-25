function [X, T, alpha] = ch3_col_unpack(z, p)
%CH3_COL_UNPACK  Split the collocation decision vector.
%
%   [X, T, alpha] = ch3_col_unpack(z, p)
%
% Inverse of ch3_col_pack.  N is inferred from the length of z, so the same
% code works for any node count without threading it through separately.
%
% See also CH3_COL_PACK.

n_alpha = p.ny * p.n_ctrl;
nX      = numel(z) - 1 - n_alpha;

if mod(nX, p.nx) ~= 0
    error('ch3_col_unpack:size', ...
          'z length %d is not consistent with nx = %d and %d alpha entries.', ...
          numel(z), p.nx, n_alpha);
end

N = nX / p.nx;

X     = reshape(z(1:nX), p.nx, N);
T     = z(nX+1);
alpha = reshape(z(nX+2:end), p.ny, p.n_ctrl);

end
