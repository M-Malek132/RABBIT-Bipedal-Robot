function [X, T, alpha, theta_pm] = ch3_col_unpack(z, p)
%CH3_COL_UNPACK  Split the collocation decision vector.
%
%   [X, T, alpha]           = ch3_col_unpack(z, p)
%   [X, T, alpha, theta_pm] = ch3_col_unpack(z, p)
%
% Inverse of ch3_col_pack.  N is inferred from the length of z, so the same
% code works for any node count without threading it through separately.
%
% theta_pm = [theta_minus; theta_plus]: the last two entries of z when
% p.free_theta is set, otherwise p's own fixed values. Code that evaluates the
% gait must use these rather than p's -- ch3_col_effective_params does that.
%
% See also CH3_COL_PACK, CH3_COL_EFFECTIVE_PARAMS.

n_alpha = p.ny * p.n_ctrl;
free    = isfield(p, 'free_theta') && ~isempty(p.free_theta) && p.free_theta;
n_extra = 2 * free;
nX      = numel(z) - 1 - n_alpha - n_extra;

if nX <= 0 || mod(nX, p.nx) ~= 0
    hint = '';
    if free && mod(numel(z) - 1 - n_alpha, p.nx) == 0
        hint = [' It has the length of a FIXED-theta vector: pass it through ' ...
                'ch3_col_theta_augment first.'];
    elseif ~free && mod(numel(z) - 3 - n_alpha, p.nx) == 0
        hint = [' It has the length of a FREE-theta vector: this p needs ' ...
                'p.free_theta = true (the p saved with the gait carries it).'];
    end
    error('ch3_col_unpack:size', ...
          'z length %d is not consistent with nx = %d and %d alpha entries%s.%s', ...
          numel(z), p.nx, n_alpha, ...
          ternary(free, ' plus [theta_minus; theta_plus]', ''), hint);
end

N = nX / p.nx;

X     = reshape(z(1:nX), p.nx, N);
T     = z(nX+1);
alpha = reshape(z(nX+2 : nX+1+n_alpha), p.ny, p.n_ctrl);

if nargout > 3
    if free
        theta_pm = z(end-1:end);
        theta_pm = theta_pm(:);
    else
        theta_pm = [p.theta_minus; p.theta_plus];
    end
end

end

% ---------------------------------------------------------------------------
function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end
