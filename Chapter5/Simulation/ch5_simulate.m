function out = ch5_simulate(p, varargin)
%CH5_SIMULATE  Run one Chapter-5 study: sampled-data control over [0, T].
%
%   out = ch5_simulate(p)
%   out = ch5_simulate(p, 'x0', x0, 'verbose', false)
%
% One loop, one control rate, no hybrid events. At each sample:
%
%   1. evaluate the controller at the measured state          ch5_control
%   2. hold u and integrate the plant over one period         ch5_ode_rhs
%
% ------------------------------------------------- the admissibility gate
% Before integrating anything, Corollary 5.2 is checked at the ACTUAL x0 for
% the ECBF controllers. This is not decoration. The exponential envelope
%
%       h(t) >= Cb exp(Ab t) eta_b(0)
%
% is only nonnegative when eta_b(0) sits in the chain of sets C_0..C_rb, and
% whether it does depends on the poles AND on x0 together. Running from an
% inadmissible x0 produces a controller that enforces its row faithfully at
% every instant and still lets h go negative -- which looks exactly like a bug
% in the barrier and is not one. p.ecbf.admissibility decides whether that is
% an error, a warning, or ignored.
%
% ---------------------------------------------------------------- what is logged
% Everything the chapter's figures need, at the CONTROL rate rather than the
% integrator's: h, eta_b, u, mu, the CLF value, the QP slack, whether the
% barrier row was active, and y_rb (Remark 5.6). The state is logged at the
% control rate too -- these plants have no fast transients between samples and
% a uniform grid makes the runs directly comparable.
%
% Inputs
%   p : parameter struct
%   name/value:
%     'x0'      override the initial state
%     'verbose' progress line (default false)
%
% Outputs
%   out : struct with
%           .t (1xN) .x (nx x N) .u (nu x N) .mu (ny x N)
%           .h (1xN) .eta_b (rb x N) .V .delta .y_rb
%           .cbf_active .feasible   logical rows
%           .h_min .h_min_t .violated
%           .adm      the admissibility report (empty for non-ECBF runs)
%           .p .ok .reason .cpu
%
% See also CH5_CONTROL, CH5_ECBF_ADMISSIBLE, CH5_MAIN, CH5_REPORT.

o = struct('x0', [], 'verbose', false);
for k = 1:2:numel(varargin), o.(lower(varargin{k})) = varargin{k+1}; end

x0 = o.x0;
if isempty(x0), x0 = ch5_x0(p); end
x0 = x0(:);

nx = p.sys.nx;  nu = p.sys.nu;  ny = p.sys.ny;
dt = p.control_dt;
N  = round(p.T / dt) + 1;

is_ecbf = startsWith(lower(p.controller), 'ecbf');

%% --------------------------------------------- Corollary 5.2, before anything
adm = [];
if is_ecbf
    b0 = ch5_barrier(x0, p);
    e  = ch5_ecbf_gain(p, b0.rb);
    adm = ch5_ecbf_admissible(x0, p, e);
    switch lower(p.ecbf.admissibility)
        case 'error'
            if ~adm.ok
                error('ch5_simulate:inadmissible', ...
                      'Corollary 5.2 fails at x0: %s', adm.msg);
            end
        case 'warn'
            if ~adm.ok
                fprintf(['  [Cor 5.2] %s\n' ...
                         '            the invariance claim is NOT supported ' ...
                         'from this x0.\n'], adm.msg);
            end
        case 'off'
            % caller has taken responsibility
        otherwise
            error('ch5_simulate:admissibilityMode', ...
                  'p.ecbf.admissibility must be error|warn|off.');
    end
else
    b0 = ch5_barrier(x0, p);
    e  = [];
end
rb = b0.rb;

%% ------------------------------------------------------------------ storage
t     = (0:N-1) * dt;
X     = nan(nx, N);
U     = nan(nu, N);
MU    = nan(ny, N);
H     = nan(1,  N);
ETAB  = nan(rb, N);
Vv    = nan(1,  N);
DEL   = nan(1,  N);
YRB   = nan(1,  N);
KBE   = nan(1,  N);
ACT   = false(1, N);
FEAS  = true(1,  N);

ok     = true;
reason = 'completed';
tic;

x = x0;
for i = 1:N

    if ~all(isfinite(x))
        ok = false;
        reason = sprintf('non-finite state at t = %.4f s', t(i));
        break;
    end

    [u, in] = ch5_control(x, p, e);

    X(:,i)    = x;
    U(:,i)    = u;
    MU(:,i)   = in.mu;
    H(i)      = in.b.h;
    ETAB(:,i) = in.b.eta_b;
    Vv(i)     = in.qp.V;
    DEL(i)    = in.qp.delta;
    ACT(i)    = in.qp.cbf_active;
    FEAS(i)   = in.qp.feasible;
    % Kb*eta_b is logged alongside y_rb because y_rb alone cannot be judged.
    % The ECBF row is  y_rb = mu_b + Kb eta_b >= 0, and on the pendulum
    % Kb eta_b runs to several thousand -- so a y_rb of -3e-5 is the solver
    % sitting on the constraint to its own tolerance, not the guarantee
    % failing. ch5_report divides by this before deciding.
    if isfield(in.qp, 'y_rb'),   YRB(i) = in.qp.y_rb;   end
    if isfield(in.qp, 'Kb_eta'), KBE(i) = in.qp.Kb_eta; end

    if i == N, break; end

    x = advance(x, u, dt, p);

    if o.verbose && mod(i, max(1, round(N/10))) == 0
        fprintf('    %5.1f%%  t = %6.2f s   h = %+8.4f\n', ...
                100*i/N, t(i), H(i));
    end
end

cpu = toc;

%% ------------------------------------------------------------------ summary
keep = ~isnan(H);
[h_min, imin] = min(H(keep));
tk = t(keep);
if isempty(h_min), h_min = NaN; h_min_t = NaN; else, h_min_t = tk(imin); end

out = struct('t', t, 'x', X, 'u', U, 'mu', MU, 'h', H, 'eta_b', ETAB, ...
             'V', Vv, 'delta', DEL, 'y_rb', YRB, 'Kb_eta', KBE, ...
             'cbf_active', ACT, 'feasible', FEAS, ...
             'h_min', h_min, 'h_min_t', h_min_t, 'violated', h_min < 0, ...
             'adm', adm, 'p', p, 'ok', ok, 'reason', reason, 'cpu', cpu, ...
             'n', sum(keep));

end

% ---------------------------------------------------------------------------
function x = advance(x, u, dt, p)
%ADVANCE  One control period with u held.
switch lower(p.integrator)

    case 'ode45'
        opts = odeset('RelTol', p.ode_opts.RelTol, 'AbsTol', p.ode_opts.AbsTol);
        [~, Z] = ode45(@(tt,zz) ch5_ode_rhs(tt, zz, u, p), [0 dt], x, opts);
        x = Z(end,:).';

    case 'rk4'
        m = p.n_substeps;
        hs = dt / m;
        for j = 1:m
            k1 = ch5_ode_rhs(0, x,           u, p);
            k2 = ch5_ode_rhs(0, x + hs/2*k1, u, p);
            k3 = ch5_ode_rhs(0, x + hs/2*k2, u, p);
            k4 = ch5_ode_rhs(0, x + hs*k3,   u, p);
            x  = x + (hs/6)*(k1 + 2*k2 + 2*k3 + k4);
        end

    otherwise
        error('ch5_simulate:integrator', ...
              'p.integrator must be ode45|rk4 (got "%s").', p.integrator);
end
end
