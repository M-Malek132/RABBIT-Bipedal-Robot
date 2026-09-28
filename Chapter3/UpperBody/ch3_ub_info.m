function [active, optimize, n_vars] = ch3_ub_info(p)
%CH3_UB_INFO  Is an upper-body wrench switched on, and is it being optimized?
%
%   [active, optimize, n_vars] = ch3_ub_info(p)
%
% The single place that reads p.ub's switches, so every patched Chapter-3 file
% (control_affine, zd_point, col_pack/unpack/bounds/cost/effective_params)
% agrees on the answer. With no p.ub field -- every gait saved before this
% feature -- all three outputs are false/0 and NOTHING in Chapter 3 changes.
%
% Outputs
%   active   : true when p.ub.mode is 'profile' or 'bezier'
%   optimize : true when the Bezier coefficients p.ub.beta are appended to the
%              collocation decision vector z (requires mode 'bezier')
%   n_vars   : number of entries appended to z (3 * n_ub_ctrl, or 0)
%
% See also CH3_UB_DEFAULTS, CH3_UB_WRENCH.

active = false;  optimize = false;  n_vars = 0;
if ~isfield(p, 'ub') || isempty(p.ub) || ~isfield(p.ub, 'mode')
    return;
end
switch lower(p.ub.mode)
    case 'none'
        return;
    case {'profile', 'bezier'}
        active = true;
    otherwise
        error('ch3_ub_info:mode', ...
              'p.ub.mode must be ''none'', ''profile'' or ''bezier'' (got ''%s'').', ...
              p.ub.mode);
end
optimize = isfield(p.ub, 'optimize') && ~isempty(p.ub.optimize) && p.ub.optimize;
if optimize
    if ~strcmpi(p.ub.mode, 'bezier')
        error('ch3_ub_info:optimizeMode', ...
              'p.ub.optimize needs p.ub.mode = ''bezier'' (a fixed profile has nothing to solve for).');
    end
    n_vars = numel(p.ub.beta);
end
end
