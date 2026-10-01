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
% WITH AN OPTIMIZED UPPER-BODY WRENCH z ENDS IN beta, AFTER theta
% (ch3_col_pack: [X(:); T; alpha(:); theta_pm; beta(:)]). The pair is inserted
% in front of that tail and the length test is made without it. Testing the
% whole length fails on every such z (1174 = 14*81 + 1 + 24 alpha + 15 beta
% fits neither layout), and appending would put theta where unpack reads beta.
%
% See also CH3_COL_PACK, CH3_COL_SOLVE.

free = isfield(p, 'free_theta') && ~isempty(p.free_theta) && p.free_theta;
if ~free
    return;
end

n_alpha = p.ny * p.n_ctrl;
[~, ~, n_ub] = ch3_ub_info(p);           % upper-body beta, stored last
z = z(:);
n_head = numel(z) - n_ub;
if mod(n_head - 1 - n_alpha, p.nx) == 0
    z = [z(1:n_head); p.theta_minus; p.theta_plus; z(n_head+1:end)];
elseif mod(n_head - 3 - n_alpha, p.nx) ~= 0
    error('ch3_col_theta_augment:size', ...
          'z length %d fits neither a fixed- nor a free-theta layout.', numel(z));
end

end
