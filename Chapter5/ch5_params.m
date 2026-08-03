function p = ch5_params(varargin)
%CH5_PARAMS  Single source of truth for the Chapter-5 (safety) pipeline.
%
%   p = ch5_params()                        defaults for the spring-mass study
%   p = ch5_params('system', 'pendulum')    defaults for the pendulum study
%   p = ch5_params('name', value, ...)      with overrides (dotted names OK)
%
% Chapter 5 is DELIBERATELY NOT ABOUT RABBIT. The chapter's own words: the
% method is validated on a relative-degree-6 linear system and a
% relative-degree-4 nonlinear system, and only then ("in the next Chapter")
% applied to 3D walking. So this chapter carries its own two plants and shares
% nothing with Chapter 3/4 but the code style.
%
% ---------------------------------------------------------------- the two ideas
% Chapter 3 asked a QP for STABILITY (a CLF row). Chapter 5 adds a second row
% asking for SAFETY -- forward invariance of C = {x : h(x) >= 0} -- and the
% whole chapter is about what that second row has to look like.
%
%   Section 5.1  RECIPROCAL CBF, B = 1/h, condition Bdot <= gamma/B.
%                Works, and only works, when h has RELATIVE DEGREE 1. Kept here
%                in full because its failure at higher relative degree is the
%                chapter's motivation, and a failure you can run is worth more
%                than one you can quote (see ch5_test_qp, check 5).
%
%   Section 5.2  EXPONENTIAL CBF. h itself is the barrier, and the condition is
%                placed on its rb-th derivative via pole placement:
%                     h^(rb) >= -Kb*eta_b,    eta_b = [h; hdot; ...; h^(rb-1)].
%                Arbitrary relative degree, and the poles are the design knob.
%
% ------------------------------------------------------------ NAMING WARNING
% The thesis writes "B(x)" for BOTH the reciprocal barrier of (5.6), where the
% safe set is {B < inf} and B blows up at the boundary, AND for the exponential
% barrier of Definition 5.1, where the safe set is {B >= 0} and B is the
% constraint itself. They are different objects with opposite sign conventions
% and mixing them is a sign bug that still satisfies every dimension check.
%
% This code therefore never uses the bare name "B". It uses:
%
%   h(x)     the constraint function.  SAFE SET IS ALWAYS  C = {h(x) >= 0}.
%   B_rec    = 1/h,  the Section 5.1 reciprocal barrier, only ever built when
%            h has relative degree 1.
%   eta_b    = [h; hdot; ...; h^(rb-1)],  eq (5.10), the Section 5.2 object.
%
% See also CH5_SYSTEM, CH5_BARRIER, CH5_ECBF_GAIN, CH5_MAIN.

p = struct();

%% ------------------------------------------------------------------ plant
% 'springmass'  Fig. 5.1, eq (5.32)-(5.35). Linear, 3 DOF, 1 input,
%               relative degree 6 from u to x3.
% 'pendulum'    Fig. 5.2, eq (5.36)-(5.37). Nonlinear, 4 DOF, 2 inputs,
%               relative degree 4 from tau to theta.
p.system = 'springmass';

%% ------------------------------------------------------------- controller
%   'clfqp'      Chapter-3 style CLF-QP, no barrier row at all. This is the
%                controller that VIOLATES the constraint in Fig. 5.3a / 5.4a,
%                and it is in the chapter to establish that the constraint is
%                not satisfied for free.
%   'cbfclfqp'   Section 5.1, eq (5.7)/(5.8). Reciprocal barrier.
%   'cbfclfqp_viol'  Section 5.1 again, but routed through the Virtual
%                Input-Output Linearization of (5.15). Remark 5.4 claims this
%                gives the IDENTICAL control; ch5_test_qp checks it.
%   'ecbfclfqp'  Section 5.2, eq (5.31). The chapter's result.
p.controller = 'clfqp';

%% ------------------------------------------- CLF (the stability row, 3.19)
p.clf = struct();

% Q of the CARE / Lyapunov equation. Sized r*ny at build time when left empty.
p.clf.Q = [];

% 'care'  F'P + PF - P G G' P + Q = 0. Pointwise feasible for every eta, so the
%         stability row never makes the QP infeasible on its own.
% 'lyap'  A'P + PA = -Q for a pole-placed A. Kept for parity with ch3_res_clf.
p.clf.construction = 'care';

% Poles used only by the 'lyap' construction (r per output).
p.clf.poles = [];

% eps of the RES-CLF eta-scaling, I_eps = blkdiag(I/eps^(r-1), ..., I).
% eps = 1 is the plain CLF, which is what (5.31) is written with.
p.clf.eps = 1.0;

% lambda of (5.31)'s  Vdot + lambda*V <= delta.  Empty means "use the rate the
% CARE actually guarantees", c3/eps with c3 = lam_min(Q)/lam_max(P) -- i.e. the
% Chapter-3 convention. Set a number to demand a specific rate.
p.clf.lambda = [];

% p of (5.31)'s penalty p*delta^2 on the stability slack.
%
% THE SLACK IS ON THE CLF ROW ONLY. The barrier row is never relaxed: a safety
% constraint you are willing to violate for a smoother input is not a safety
% constraint. That asymmetry is the entire reason the chapter's QP has both
% rows rather than one weighted objective.
p.clf.slack_penalty = 1e6;

%% ---------------------------------------------- CBF, Section 5.1 (5.3)/(5.6)
p.cbf = struct();

% gamma of  Bdot <= gamma/B  with B = 1/h, eq (5.3) with the reciprocal
% candidate (5.6). Equivalent to  hdot >= -gamma*h^3  (see ch5_barrier).
p.cbf.gamma = 1.0;

%% ----------------------------------------- Exponential CBF, Section 5.2
p.ecbf = struct();

% DESIRED POLES p_b of (5.28). Must be rb of them, all real and STRICTLY
% POSITIVE -- they enter as -p_i, and Theorem 5.2 needs A_b Hurwitz AND "total
% negative" (real negative eigenvalues, not merely negative real parts). A
% complex pair would satisfy Hurwitz and break Proposition 5.1, because the
% recursion (5.30) is a chain of REAL first-order filters.
%
% Empty means "use the per-system default in ch5_system".
p.ecbf.poles = [];

% Corollary 5.2 says the poles are not free: p_i >= -ydot_{i-1}(x0)/y_{i-1}(x0).
% ch5_ecbf_admissible checks it at the actual initial condition.
%   'error'  refuse to run from an x0 the poles do not cover
%   'warn'   report and continue (the invariance claim is then unsupported)
%   'off'    do not check
p.ecbf.admissibility = 'warn';

%% ------------------------------------------------- the constraint to enforce
% Override for the plant's own default safety level, which lives in
% p.plant.constraint once ch5_system has run. This is how the figure sweeps
% read: Fig. 5.3 varies x3max, Fig. 5.4 varies p2min, everything else fixed.
%
%   []              use the plant default
%   a scalar        replace only the LEVEL (x3max, or p2min)
%   a struct        replace the whole spec, .type and .value
%
% Constraint types understood by ch5_barrier:
%   springmass  'x3_max'  x3 <= value        relative degree 6
%               'v1_max'  xdot1 <= value     relative degree 1  <- see below
%   pendulum    'py_min'  py2 >= value       relative degree 4
%
% 'v1_max' is not in the thesis. It is here because Section 5.1's reciprocal
% CBF is only defined for relative degree 1, so without a relative-degree-1
% constraint on one of these plants there is nowhere to demonstrate that
% Section 5.1 WORKS -- only places it fails. Remark 5.4's claim that the VIOL
% form (5.15) and the direct form (5.7) give identical controls is checked on
% it too, since that claim is likewise vacuous where the barrier row is
% uncontrollable.
p.constraint = [];

%% ------------------------------------------------------------ input limits
% Ac(x) u <= bc(x) of (5.8)/(5.31), in the only form either plant needs: a
% symmetric box. Empty disables it, which is the chapter's own setting for both
% validation systems -- Fig. 5.4c in particular needs several hundred Nm and
% would be a different experiment with a box on it.
p.u_max = [];

%% ------------------------------------------------------------- integration
% Sampled-data control, for the same reason ch4_params gives: a QP solution is
% only piecewise smooth in x, and an adaptive ODE solver asked to integrate
% through that will shrink its step forever chasing a discontinuity in the
% derivative that is really a change of active set.
p.control_dt = 5e-3;

% Per-control-interval integrator. 'ode45' is adaptive and accurate; 'rk4' is
% fixed-step with p.n_substeps substeps and is ~8x faster, used by the tests.
p.integrator = 'ode45';
p.n_substeps = 10;
p.ode_opts   = struct('RelTol', 1e-9, 'AbsTol', 1e-11);

% Total simulated time [s]. Overridden per system in ch5_system.
p.T = [];

%% -------------------------------------------------------------- overrides
for k = 1:2:numel(varargin)
    p = set_field(p, varargin{k}, varargin{k+1});
end

%% ------------------------------------------------- pull in the plant block
% Done AFTER the overrides so that p.system can be overridden, but the plant
% block it selects is then merged WITHOUT clobbering anything the caller set
% explicitly -- so ch5_params('system','pendulum','T',5) keeps T = 5.
explicit = varargin(1:2:end);
p = ch5_system(p, explicit);

end

% ---------------------------------------------------------------------------
function s = set_field(s, name, value)
parts = strsplit(name, '.');
if ~isfield(s, parts{1})
    error('ch5_params:unknownField', 'Unknown parameter "%s".', name);
end
if numel(parts) == 1
    s.(parts{1}) = value;
else
    s.(parts{1}) = set_field(s.(parts{1}), strjoin(parts(2:end), '.'), value);
end
end
