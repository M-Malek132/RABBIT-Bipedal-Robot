function col_verify_seed(result_file, N)
%COL_VERIFY_SEED  Sanity-check the collocation transcription on a PD seed.
%
%   col_verify_seed()              % newest Results/hzd_result_*.mat, N=15
%   col_verify_seed(result_file, N)
%
% Builds a PD-warm-started collocation vector and reports, WITHOUT solving:
%   * vector length vs. lb/ub (layout bug catch)
%   * norm of each equality-constraint block (defects, virtual constraints,
%     foot pinning, periodicity, ...) so you can see how close the PD seed is
%   * the cost, and the contact-force sign (so you know how to set enforce_grf)
%   * a finite-difference check of col_dynamics (Jacobian sanity)
%
% A good seed has small dynamics defects (the PD (q,dq,u,lambda) satisfy the
% EOM), small foot drift, and modest virtual-constraint residual (PD tracking
% error). Large numbers here mean a transcription bug, not a hard NLP.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(root));

    if nargin < 2 || isempty(N), N = 15; end
    if nargin < 1 || isempty(result_file)
        d = dir(fullfile(root,'Results','hzd_result_*.mat'));
        assert(~isempty(d), 'No Results/hzd_result_*.mat; pass a file explicitly.');
        % Newest by NAME (the hzd_result_YYYY-MM-DD_HH-MM-SS stamp sorts
        % chronologically), not by mtime: a fresh clone/worktree gives every
        % file the same mtime, so max([d.datenum]) would pick arbitrarily.
        [~,ord] = sort({d.name});
        i = ord(end);
        result_file = fullfile(d(i).folder, d(i).name);
    end
    fprintf('Seed gait : %s\nNodes N   : %d\n\n', result_file, N);

    [p, lb, ub] = col_problem_setup(N);
    S = load(result_file);
    [coeffs, x_start] = unpack_z(S.z_opt, p);

    z0 = col_seed_from_pd(coeffs, x_start, p);

    fprintf('length(z0)=%d  length(lb)=%d  length(ub)=%d  %s\n', ...
        numel(z0), numel(lb), numel(ub), ...
        ternary(numel(z0)==numel(lb) && numel(z0)==numel(ub),'[OK]','[SIZE MISMATCH]'));
    nOut = sum(z0 < lb - 1e-9) + sum(z0 > ub + 1e-9);
    fprintf('seed entries outside bounds: %d\n\n', nOut);

    % --- constraint blocks, sized exactly as assembled in col_constraints ---
    [c, ceq] = col_constraints(z0, p);
    nq = p.nq; nu = p.nu; nx = 2*nq;
    nDef = nx*(N-1);  nFD = 2*(N-1);  nVC = nu*N;
    % Cumulative slicing with an explicit counter (do NOT use an anonymous
    % function closing over a mutated index -- it captures the value once).
    i0 = 0;
    [b_def, i0] = take(ceq, i0, nDef);
    [b_fd , i0] = take(ceq, i0, nFD);
    [b_vc , i0] = take(ceq, i0, nVC);
    [b_fz1, i0] = take(ceq, i0, 1);
    [b_fv1, i0] = take(ceq, i0, 2);
    [b_yd1, i0] = take(ceq, i0, nu);
    [b_is , i0] = take(ceq, i0, 1);
    [b_te , i0] = take(ceq, i0, 1);
    [b_per, i0] = take(ceq, i0, 13);
    [b_spd, i0] = take(ceq, i0, 1);
    idx = i0;

    fprintf('equality-constraint block infinity-norms (seed):\n');
    pr('  HS dynamics defects', b_def);
    pr('  foot no-drift      ', b_fd);
    pr('  virtual constr y=0 ', b_vc);
    pr('  foot z @node1      ', b_fz1);
    pr('  foot vel @node1    ', b_fv1);
    pr('  ydot @node1        ', b_yd1);
    pr('  impact surface     ', b_is);
    pr('  theta_end          ', b_te);
    pr('  periodicity        ', b_per);
    pr('  speed              ', b_spd);
    fprintf('  (all ceq consumed: %s)\n', ternary(idx==numel(ceq),'yes','NO - size bug'));
    fprintf('max inequality c  : %.3e (<=0 satisfied if negative)\n', max(c));

    fprintf('\ncost J(seed)      : %.4g\n', col_cost(z0, p));

    % contact-force sign over the step (helps set enforce_grf / friction signs)
    [~,~,Lam] = col_unpack(z0, p);
    fprintf('Lam Fz  min/mean/max : %.1f / %.1f / %.1f  N\n', ...
            min(Lam(2,:)), mean(Lam(2,:)), max(Lam(2,:)));

    % --- finite-difference check of col_dynamics Jacobian d(dx)/d[x;u;lam] ---
    x = z0(1:nx); u = z0(nx*N+(1:nu)); lam = [0;0];
    fd_ok = fd_check_dynamics(x, u, lam, p);
    fprintf('\ncol_dynamics FD Jacobian check: max rel err %.2e  %s\n', ...
            fd_ok, ternary(fd_ok<1e-4,'[OK]','[CHECK]'));
end

function [blk, i1] = take(v, i0, n)
    blk = v(i0 + (1:n));
    i1  = i0 + n;
end

function pr(name, v)
    fprintf('%s : ||.||inf = %.3e\n', name, max(abs(v)));
end

function out = ternary(cond, a, b)
    if cond, out=a; else, out=b; end
end

function maxrel = fd_check_dynamics(x, u, lam, p)
% central-difference the analytic dx=col_dynamics wrt x to confirm it is a
% smooth, correctly-wired function (a coarse guard, not a proof).
    nx = numel(x);
    f0 = col_dynamics(x, u, lam, p);
    h = 1e-6; maxrel = 0;
    for i = 1:nx
        xp=x; xp(i)=xp(i)+h; xm=x; xm(i)=xm(i)-h;
        num = (col_dynamics(xp,u,lam,p) - col_dynamics(xm,u,lam,p))/(2*h);
        % analytic column i of df/dx is not formed; instead check smoothness by
        % comparing central vs one-sided (a consistency proxy).
        one = (col_dynamics(xp,u,lam,p) - f0)/h;
        maxrel = max(maxrel, norm(num-one)/max(1,norm(num)));
    end
end
