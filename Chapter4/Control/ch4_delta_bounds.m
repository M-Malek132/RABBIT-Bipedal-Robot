function B = ch4_delta_bounds(X, alpha, p, scales, opts)
%CH4_DELTA_BOUNDS  Measure Delta1, Delta2 so the bounds (4.10) can be CHOSEN.
%
%   B = ch4_delta_bounds(X, alpha, p)
%   B = ch4_delta_bounds(X, alpha, p, scales)
%   B = ch4_delta_bounds(X, alpha, p, scales, opts)
%
% The robust CLF-QP needs Delta1max and Delta2max in (4.10), and its entire
% guarantee is conditional on those numbers actually bounding the uncertainty.
% Picking them by eye is therefore not a shortcut, it is the difference between
% a guarantee and a decoration. This function measures the induced uncertainty
% (4.4) along a trajectory, over whatever set of model perturbations you intend
% to survive, and reports what the bounds have to be.
%
% WHY MEASURING ON THE ORBIT ALONE IS NOT ENOUGH.  The controller does not run
% on the nominal orbit -- if it did there would be no tracking error to correct.
% It runs in a neighborhood, and Delta varies over that neighborhood. So the
% sweep optionally jitters each sampled state (opts.jitter, opts.n_jitter) and
% takes the worst case over the cloud. A bound fitted only to the orbit is
% systematically optimistic, and the failure mode it produces -- an infeasible
% or wrongly-certified QP off the orbit -- is exactly where you need the bound
% to hold.
%
% READING THE RESULT.  Delta2 is the one to watch. ch4_ctrl_rclf_qp is
% pointwise feasible only while Delta2max < 1 (see its header): at Delta2max = 1
% the worst-case model can null the control's effect on V entirely, and beyond
% it the model can reverse it. So a measured n2 near or above 1 is not a number
% to pad with a safety factor, it is the method reporting that this much
% uncertainty is outside what a worst-case bound can cover.
%
% Inputs
%   X      : 14 x nt state trajectory to sample along (e.g. collocation nodes,
%            or a simulated rollout -- denser is better)
%   alpha  : ny x n_ctrl virtual constraint coefficients
%   p      : parameter struct
%   scales : vector of mass_scale values to cover. Default the chapter's
%            Cases I-III, [1 0.7 1.5]; add 3 for Case IV.
%   opts   : struct with
%              .jitter    std-dev of the state perturbation (default 0.05)
%              .n_jitter  cloud size per sample (default 0, i.e. orbit only)
%              .safety    multiplier applied to the measured max (default 1.25)
%              .loads     vector of load_mass values to also cover (default 0)
%              .verbose   print the per-scale table (default true)
%
% Output
%   B : struct with
%         .delta1_max .delta2_max   the RECOMMENDED bounds (measured * safety)
%         .n1_max .n2_max           raw measured maxima over everything
%         .per_scale                struct array, one per (scale, load) pair,
%                                   with .mass_scale .load_mass .n1 .n2 .n2_scalar
%         .feasible                 false if delta2_max >= 1
%
% See also CH4_UNCERTAINTY, CH4_CTRL_RCLF_QP.

if nargin < 4 || isempty(scales), scales = [1 0.7 1.5]; end
if nargin < 5, opts = struct(); end

opts = fill_defaults(opts, struct('jitter', 0.05, 'n_jitter', 0, ...
                                  'safety', 1.25, 'loads', 0, ...
                                  'verbose', true));

nt = size(X, 2);
rng(7);                      % reproducible jitter cloud

combos = struct('mass_scale', {}, 'load_mass', {}, ...
                'n1', {}, 'n2', {}, 'n2_scalar', {});

if opts.verbose
    fprintf('\n%s\n MEASURED MODEL UNCERTAINTY (eq 4.4) over %d states', ...
            repmat('=',1,66), nt);
    if opts.n_jitter > 0
        fprintf(' x %d jitters (sd %.3f)', opts.n_jitter, opts.jitter);
    end
    fprintf('\n%s\n', repmat('-',1,66));
    fprintf(' %10s %10s %12s %10s %10s\n', ...
            'mass sc.', 'load kg', '|Delta1|', '|Delta2|', 'iso(D2)');
    fprintf('%s\n', repmat('-',1,66));
end

for is = 1:numel(scales)
    for il = 1:numel(opts.loads)

        unc = struct('mass_scale', scales(is), 'load_mass', opts.loads(il));

        n1 = 0; n2 = 0; n2s = 0;

        for k = 1:nt
            % the sampled state itself, plus its jitter cloud
            for j = 0:opts.n_jitter
                xk = X(:, k);
                if j > 0
                    xk = xk + opts.jitter * randn(size(xk));
                end

                D = ch4_uncertainty(xk, alpha, p, unc);

                % A near-singular NOMINAL decoupling matrix makes Delta2
                % arbitrarily large for reasons that have nothing to do with
                % model error -- it is the inversion blowing up, not the
                % uncertainty growing. Skip those rather than let one bad pose
                % set a bound the controller then has to carry everywhere.
                if ~isfinite(D.n1) || ~isfinite(D.n2), continue; end

                n1  = max(n1,  D.n1);
                n2  = max(n2,  D.n2);
                n2s = max(n2s, D.n2_scalar);
            end
        end

        combos(end+1) = struct('mass_scale', scales(is), ...
                               'load_mass',  opts.loads(il), ...
                               'n1', n1, 'n2', n2, 'n2_scalar', n2s); %#ok<AGROW>

        if opts.verbose
            fprintf(' %10.2f %10.1f %12.2f %10.4f %10.4f\n', ...
                    scales(is), opts.loads(il), n1, n2, n2s);
        end
    end
end

B = struct();
B.per_scale  = combos;
B.n1_max     = max([combos.n1]);
B.n2_max     = max([combos.n2]);
B.delta1_max = opts.safety * B.n1_max;
B.delta2_max = opts.safety * B.n2_max;
B.feasible   = B.delta2_max < 1;

if opts.verbose
    fprintf('%s\n', repmat('-',1,66));
    fprintf(' recommended (x%.2f safety):  delta1_max = %.1f,  delta2_max = %.4f\n', ...
            opts.safety, B.delta1_max, B.delta2_max);
    if ~B.feasible
        fprintf([' *** delta2_max >= 1: the worst-case model inside this bound can\n' ...
                 ' *** cancel the control authority over V, so the robust CLF-QP is\n' ...
                 ' *** NOT pointwise feasible. Shrink the uncertainty set or accept\n' ...
                 ' *** that this operating point is outside worst-case robustness.\n']);
    end
    fprintf('%s\n\n', repmat('=',1,66));
end

end

% ---------------------------------------------------------------------------
function s = fill_defaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i}))
        s.(f{i}) = d.(f{i});
    end
end
end
