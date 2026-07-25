function [x_all, t_all] = animate_col_result(col_file, nSteps, opts)
%ANIMATE_COL_RESULT  Animate a direct-collocation gait (Results/col_result_*.mat).
%
%   animate_col_result()                          % newest result, 4 steps
%   animate_col_result(col_file, nSteps, opts)
%   [x_all, t_all] = animate_col_result(...)
%
% A collocation solution is NOT a time series: it is N node states plus a step
% time T, with the dynamics enforced as Hermite-Simpson defects BETWEEN the
% nodes (col_constraints). Playing the N=20 nodes back directly would give a
% 20-frame step that visibly snaps between poses, so this reconstructs the
% trajectory the transcription actually implies -- on each interval the cubic
% Hermite polynomial through (x_k, f_k) and (x_{k+1}, f_{k+1}), which is exactly
% the curve whose midpoint slope the Simpson defect pins to the dynamics -- and
% samples it at opts.n_sub points per interval. No ODE solve is involved, so
% what you see is the optimizer's own answer, not a re-simulation of it.
%
% MULTI-STEP. The solution is one step of a PERIODIC orbit: col_constraints
% enforces reset(impact(x_N)) = x_1 on coordinates 2:14, with only the
% horizontal base coordinate px free to advance. Step k is therefore step 1
% translated forward by dx = px_N - px_1 (= the step length), and tiling it
% accumulates no error. The printed periodicity defect is how exactly that
% constraint was actually met -- a large value means the tiling is papering over
% a solve that did not converge, and the walk you are watching is not one the
% robot would repeat. For the CLOSED-LOOP behaviour (which for this seed orbit
% diverges within a few steps) use col_crosscheck / animate_hzd_result instead.
%
% Blue is always the STANCE leg and red the SWING leg, so at each impact the
% colours swap physical legs. Same convention as animate_hzd_result, whose
% per-step rabbit_reset_map does the relabeling explicitly.
%
% Args (all optional):
%   col_file  Results/col_result_*.mat holding z_opt,p (default: newest by name)
%   nSteps    how many periods to tile (default 4)
%   opts      struct, any subset of:
%               .n_sub    samples per node interval          (default 10)
%               .gif_file output GIF   (default Results/collocation_animation.gif)
%               .skip     draw every skip-th frame           (default 1)
%               .delay    GIF frame delay [s]                (default 0.03)
%               .title    axes title
%
% Returns the stitched trajectory (14 x nFrames) and its time vector, so you can
% feed it to other tools without re-animating.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(repo_genpath(root));

    if nargin < 2 || isempty(nSteps), nSteps = 4;        end
    if nargin < 3 || isempty(opts),   opts   = struct(); end

    % ---- resolve the result file ----------------------------------------
    % Newest by NAME, not mtime: the col_result_YYYY-MM-DD_HH-MM-SS stamp
    % sorts chronologically, while a fresh clone/worktree stamps every file
    % with the same mtime (see rabbit_hzd_collocation for the full rationale).
    if nargin < 1 || isempty(col_file)
        d = dir(fullfile(root,'Results','col_result_*.mat'));
        assert(~isempty(d), 'No Results/col_result_*.mat; run rabbit_hzd_collocation first.');
        [~,ord] = sort({d.name});
        col_file = fullfile(d(ord(end)).folder, d(ord(end)).name);
    elseif ~isfile(col_file)
        col_file = fullfile(root,'Results', col_file);
    end

    S = load(col_file);
    p = S.p;
    [X, U, Lam, ~, T] = col_unpack(S.z_opt, p);
    N  = p.N;   nx = 2*p.nq;
    h  = T/(N-1);

    % ---- node derivatives, for the Hermite slopes ------------------------
    F = zeros(nx, N);
    for k = 1:N
        F(:,k) = col_dynamics(X(:,k), U(:,k), Lam(:,k), p);
    end

    % ---- dense ONE-step trajectory from the Hermite interpolant ----------
    % tau is half-open on each interval: node k is emitted, node k+1 is left to
    % the next interval, so no pose is drawn twice. The closing node is
    % appended once at the end.
    %
    % n_sub defaults to REAL-TIME playback rather than a fixed number: the GIF
    % shows one frame every opts.delay seconds, so emitting T/delay frames per
    % step makes wall-clock playback match the gait's own step time. Asking for
    % more sub-samples than that gives a smoother but slow-motion GIF -- the
    % printed playback factor says which you got.
    delay = getfielddef(opts, 'delay', 0.03);
    m     = getfielddef(opts, 'n_sub', max(1, round(T / (delay*(N-1)))));
    m     = max(1, round(m));
    n1  = (N-1)*m;                      % frames per step, closing node excluded
    Xs  = zeros(nx, n1+1);
    ts  = zeros(1,  n1+1);
    c   = 0;
    for k = 1:N-1
        for j = 1:m
            s  = (j-1)/m;
            b0 =  2*s^3 - 3*s^2 + 1;    % x_k
            b1 =    s^3 - 2*s^2 + s;    % h*f_k
            b2 = -2*s^3 + 3*s^2;        % x_{k+1}
            b3 =    s^3 -   s^2;        % h*f_{k+1}
            c  = c + 1;
            Xs(:,c) = b0*X(:,k) + b2*X(:,k+1) + h*(b1*F(:,k) + b3*F(:,k+1));
            ts(c)   = (k-1)*h + s*h;
        end
    end
    Xs(:,end) = X(:,N);
    ts(end)   = T;

    % ---- tile nSteps periods, advancing px by one step length each -------
    dx    = X(1,N) - X(1,1);
    ntot  = nSteps*n1 + 1;
    x_all = zeros(nx, ntot);
    t_all = zeros(1,  ntot);
    for s = 1:nSteps
        idx = (s-1)*n1 + (1:n1);
        x_all(:,idx)  = Xs(:,1:n1);
        x_all(1,idx)  = x_all(1,idx) + (s-1)*dx;
        t_all(idx)    = ts(1:n1)     + (s-1)*T;
    end
    x_all(:,end) = Xs(:,end);
    x_all(1,end) = x_all(1,end) + (nSteps-1)*dx;
    t_all(end)   = ts(end)      + (nSteps-1)*T;

    % ---- report what is being tiled --------------------------------------
    xplus  = rabbit_reset_map(rabbit_impact_map(X(:,N)));
    defect = max(abs(xplus(2:14) - X(2:14,1)));
    [~, fn, ex] = fileparts(col_file);
    fprintf('collocation gait : %s%s\n', fn, ex);
    fprintf('  N = %d nodes, T = %.4f s, step = %.4f m, speed = %.4f m/s\n', ...
            N, T, dx, dx/T);
    fprintf('  %d steps, %d frames (%d sub-samples per node interval)\n', ...
            nSteps, ntot, m);
    skip     = max(1, round(getfielddef(opts,'skip',1)));
    playback = (numel(1:skip:ntot) * delay) / t_all(end);
    fprintf('  playback %.2fx real time (%.3g s delay, every %d%s frame)\n', ...
            1/playback, delay, skip, ordsuffix(skip));
    fprintf('  periodicity defect (inf-norm, coords 2:14) : %.3e\n', defect);
    if defect > 1e-3
        warning('animate_col_result:periodicity', ...
            ['periodicity defect %.2e is large -- the tiled steps are not a ' ...
             'true periodic orbit, so this animation smooths over a solve ' ...
             'that did not converge.'], defect);
    end

    % ---- hand off to the shared renderer ---------------------------------
    a.gif_file = getfielddef(opts,'gif_file', ...
                     fullfile(root,'Results','collocation_animation.gif'));
    a.skip     = skip;
    a.delay    = delay;
    a.title    = getfielddef(opts,'title', sprintf( ...
                     'RABBIT collocation gait  (%d nodes, T = %.3f s, %.2f m/s)', ...
                     N, T, dx/T));
    animate_rabbit(x_all, a);
end

% -------------------------------------------------------------------------
function v = getfielddef(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% -------------------------------------------------------------------------
function s = ordsuffix(n)
% 'st'/'nd'/'rd'/'th' for the frame-stride message.
switch mod(n,10)
    case 1,    s = 'st';
    case 2,    s = 'nd';
    case 3,    s = 'rd';
    otherwise, s = 'th';
end
if mod(n,100) >= 11 && mod(n,100) <= 13, s = 'th'; end
end

% -------------------------------------------------------------------------
function pth = repo_genpath(root)
%REPO_GENPATH  genpath(root) with dot-directories removed.
% genpath DOES descend into '.'-prefixed folders, so on a checkout that keeps
% git worktrees under .claude/ the plain genpath(root) puts a SECOND, older copy
% of every repo function on the path -- and because '.claude' sorts before
% 'Trajectory_Optimization'/'Visualization', those stale copies WIN every lookup
% (check with `which col_unpack -all`). addpath'ing this filtered list prepends
% the real working tree so it takes precedence.
    parts = strsplit(genpath(root), pathsep);
    keep  = true(size(parts));
    for i = 1:numel(parts)
        if isempty(parts{i})
            keep(i) = false;
        else
            % Test only the part BELOW root, so a repo that itself lives under a
            % dot-directory is not filtered away entirely.
            rel = parts{i}(min(numel(root)+1, numel(parts{i})+1):end);
            keep(i) = ~contains(rel, [filesep '.']);
        end
    end
    pth = strjoin(parts(keep), pathsep);
end
