function [X, PX] = ch5_extra_barriers(x, p, e)
%CH5_EXTRA_BARRIERS  The barriers p.ecbf.extra asks for, alongside p.constraint.
%
%   X        = ch5_extra_barriers(x, p)
%   X        = ch5_extra_barriers(x, p, e)    gains from e.extra, if present
%   [X, PX]  = ch5_extra_barriers(...)        and the per-barrier parameters
%
% Each entry of p.ecbf.extra is an ordinary constraint spec, so each extra
% barrier is evaluated by ch5_barrier, its gain built by ch5_ecbf_gain and its
% initial condition checked by ch5_ecbf_admissible, exactly as the main one --
% through a copy of p whose .constraint is that spec and whose .ecbf.poles are
% the spec's own (or the plant's default when the spec leaves them empty).
%
% Output
%   X  : 1 x n struct array .b (ch5_barrier output) .e (ch5_ecbf_gain output),
%        [] when p.ecbf.extra is empty
%   PX : 1 x n cell of the per-barrier parameter structs
%
% See also CH5_BARRIER, CH5_CTRL_ECBF_CLF_QP, CH5_SIMULATE.

X = []; PX = {};
if ~isfield(p, 'ecbf') || ~isfield(p.ecbf, 'extra') || isempty(p.ecbf.extra)
    return;
end
n  = numel(p.ecbf.extra);
X  = struct('b', cell(1, n), 'e', cell(1, n));
PX = cell(1, n);
for k = 1:n
    px = p;
    px.constraint = p.ecbf.extra(k);
    px.ecbf.extra = [];
    if isfield(p.ecbf.extra(k), 'poles') && ~isempty(p.ecbf.extra(k).poles)
        px.ecbf.poles = p.ecbf.extra(k).poles;
    end
    PX{k} = px;
    X(k).b = ch5_barrier(x, px);
    if nargin >= 3 && ~isempty(e) && isfield(e, 'extra') && numel(e.extra) >= k
        X(k).e = e.extra(k);
    else
        X(k).e = ch5_ecbf_gain(px, X(k).b.rb);
    end
end
end
