function z = ch3_col_pack(X, T, alpha, p, theta_pm)
%CH3_COL_PACK  Assemble the collocation decision vector.
%
%   z = ch3_col_pack(X, T, alpha, p)
%   z = ch3_col_pack(X, T, alpha, p, theta_pm)
%
% Layout:
%
%       z = [ X(:) ; T ; alpha(:) ]                      p.free_theta false
%       z = [ X(:) ; T ; alpha(:) ; theta_pm ]           p.free_theta true
%           \_ 14N _/   1   \_ ny*n_ctrl _/  \_ 2 _/
%
%   X        : 14 x N  state at every collocation node
%   T        : scalar  step duration (free -- the gait finds its own timing)
%   alpha    : ny x n_ctrl  the Bezier coefficients, which are what stage 3 is
%              actually solving for; the node states are there to make the
%              dynamics algebraic rather than an ODE inside the optimizer.
%   theta_pm : 2x1 [theta_minus; theta_plus], the phase endpoints, appended
%              only when p.free_theta is set (see ch3_params for why). Omitted,
%              it defaults to p's own values, so a fixed-theta gait packs into a
%              free-theta vector at the phase it was solved with.
%
% See also CH3_COL_UNPACK, CH3_COL_EFFECTIVE_PARAMS.

z = [X(:); T; alpha(:)];

free = isfield(p, 'free_theta') && ~isempty(p.free_theta) && p.free_theta;
if free
    if nargin < 5 || isempty(theta_pm)
        theta_pm = [p.theta_minus; p.theta_plus];
    end
    z = [z; theta_pm(:)];
end

% Upper-body wrench coefficients, LAST, when p.ub.optimize is set
% (Chapter3/UpperBody). Seeded from p.ub.beta.
[~, ub_opt, n_ub] = ch3_ub_info(p);
if ub_opt
    z = [z; p.ub.beta(:)];
end

if numel(z) ~= p.nx*size(X,2) + 1 + p.ny*p.n_ctrl + 2*free + n_ub
    error('ch3_col_pack:size', 'Packed vector has unexpected length.');
end

end
