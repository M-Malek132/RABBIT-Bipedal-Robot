function [rho, eigvals, x_start, coeffs] = col_poincare(col_file, controller)
%COL_POINCARE  Closed-loop orbital stability of a collocation gait.
%
%   [rho, eigvals] = col_poincare(col_file, controller)
%
% Extracts the spline coeffs and the node-1 state from a saved collocation
% result and computes the full-order step-to-step (Poincare) spectral radius
% under the chosen controller via the existing poincare_stability (which
% finite-differences reset(impact(simulate_hzd_gait(...))), and simulate_hzd_gait
% dispatches on p.controller). rho < 1 => the gait is orbitally STABLE.
%
% This is what distinguishes a merely-periodic collocation orbit from a
% WALKABLE one: direct collocation makes the orbit an open-loop fixed point,
% but not necessarily an attracting one under a given feedback law.
%
%   controller : 'pd' (default; ~26 fast sims) or 'clfqp' (each sim solves a
%                QP per ODE step -- SLOW, ~26 minutes-long sims; avoid).

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(root));

    if nargin < 2 || isempty(controller), controller = 'pd'; end
    if nargin < 1 || isempty(col_file)
        d = dir(fullfile(root,'Results','col_result_*.mat'));
        assert(~isempty(d), 'No Results/col_result_*.mat; run rabbit_hzd_collocation first.');
        % Newest by NAME (the col_result_YYYY-MM-DD_HH-MM-SS stamp sorts
        % chronologically), not by mtime: a fresh clone/worktree gives every
        % file the same mtime, so max([d.datenum]) would pick arbitrarily.
        [~,ord] = sort({d.name});
        i = ord(end);
        col_file = fullfile(d(i).folder, d(i).name);
    elseif ~isfile(col_file)
        col_file = fullfile(root,'Results', col_file);
    end

    S = load(col_file);
    p = S.p;
    [X, ~, ~, coeffs, ~] = col_unpack(S.z_opt, p);
    x_start = X(:,1);
    p.controller = controller;

    [rho, eigvals] = poincare_stability(coeffs, x_start, p);

    fprintf('%s\n  controller=%s   rho=%.4f   %s\n', ...
        col_file, controller, rho, ...
        ternary(rho < 1, '[STABLE]', '[UNSTABLE]'));
end

function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
