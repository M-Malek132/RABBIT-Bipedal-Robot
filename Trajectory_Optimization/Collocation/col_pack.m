function z = col_pack(X, U, Lam, coeffs, T)
%COL_PACK  Assemble the collocation decision vector from its blocks.
%   z = col_pack(X, U, Lam, coeffs, T)   (inverse of col_unpack).
    z = [X(:); U(:); Lam(:); coeffs(:); T];
end
