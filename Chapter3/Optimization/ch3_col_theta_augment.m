function z = ch3_col_theta_augment(z, p)
%CH3_COL_THETA_AUGMENT  Give a fixed-theta vector the two free-theta entries.
%
%   z = ch3_col_theta_augment(z, p)
%
% Every gait solved so far has theta_minus and theta_plus fixed, so its z has
% no room for them. To warm-start a p.free_theta solve from such a gait, append
% [p.theta_minus; p.theta_plus] -- the values it was solved at, so the starting
% point is the same gait exactly. A vector that already carries them, or any
% vector when p.free_theta is off, is returned unchanged.
%
% The two lengths cannot be confused: they differ by 2, and nx = 14 divides
% only one of them once T and alpha are taken off.
%
% See also CH3_COL_PACK, CH3_COL_SOLVE.

free = isfield(p, 'free_theta') && ~isempty(p.free_theta) && p.free_theta;
if ~free
    return;
end

n_alpha = p.ny * p.n_ctrl;
z = z(:);
if mod(numel(z) - 1 - n_alpha, p.nx) == 0
    z = [z; p.theta_minus; p.theta_plus];
elseif mod(numel(z) - 3 - n_alpha, p.nx) ~= 0
    error('ch3_col_theta_augment:size', ...
          'z length %d fits neither a fixed- nor a free-theta layout.', numel(z));
end

end
