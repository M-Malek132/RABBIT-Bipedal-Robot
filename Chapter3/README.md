# Chapter 3 — Hybrid Walking Pipeline

A from-scratch implementation of the Chapter 3 pipeline for the RABBIT planar
five-link biped:

> hybrid model → virtual constraints with a phase variable → offline
> optimization for α → I/O linearization → RES-CLF → online QP that enforces
> stability *and* physical limits

Everything here is new code. The only things reused from the surrounding repo
are the **robot model itself** — the symbolically generated `M, V, G, J_st,
J_sw, Jdotdq_*, P_st, P_sw, T1–T4, Tt` in `Dynamics/` — and the **B-spline
evaluator** (`Trajectory_Optimization/BSpline.m`), which is kept as a
cross-check on the Bézier basis.

---

## Quick start

```matlab
startup                 % from the repo root; adds Chapter3/ to the path
ch3_test_all            % verify every stage (a few minutes)
out = ch3_main;         % solve a gait end to end and report on it
```

To look at a solved gait:

```matlab
ch3_plot_gait(out.z_opt, out.p, 'Results/ch3_gait.png');
ch3_compare_controllers(out.z_opt, out.p);
```

> **Long solves and this MATLAB install.** Solves here have been observed to
> die mid-run from crashes inside MATLAB's own add-on registry and worker
> threads — unrelated to this code, but fatal to the process. `ch3_col_solve`
> therefore checkpoints `z` every `p.checkpoint_every` iterations to
> `p.checkpoint_file`. If a run dies, load the checkpoint and carry on; nothing
> is lost but the last few iterations. Running `matlab -nojvm` avoids one of
> the two observed crash paths.

---

## The eight stages, and where each one lives

| # | Stage | Files |
|---|-------|-------|
| 1 | **Hybrid model** `ẋ = f + gu`, `x⁺ = Δ(x⁻)` | `Model/ch3_control_affine.m`, `ch3_guard.m`, `ch3_impact.m`, `ch3_relabel.m` |
| 2 | **Virtual constraints** `y = y₀(q) − y_d(s(q),α)` | `VirtualConstraints/ch3_phase.m`, `ch3_bezier.m`, `ch3_yd.m`, `ch3_outputs.m` |
| 3 | **Offline optimization** for α | `Optimization/ch3_col_*.m`, `ch3_seed.m` |
| 4 | **I/O linearization** `ÿ = L_f²y + L_gL_fy·u` | `Control/ch3_io_lin.m` |
| 5 | **PD baseline** | `Control/ch3_ctrl_pd.m` |
| 6 | **RES-CLF** | `Control/ch3_res_clf.m`, `ch3_clf_eval.m` |
| 7 | **CLF-QP** | `Control/ch3_ctrl_clf_qp.m` (`constrained = false`) |
| 8 | **Constrained CLF-QP** | `Control/ch3_ctrl_clf_qp.m` (`constrained = true`) |

Supporting: `Simulation/` (`ch3_ode_rhs`, `ch3_step`, `ch3_simulate`),
`Analysis/` (`ch3_report`, `ch3_forces`, `ch3_poincare`, `ch3_plot_gait`,
`ch3_compare_controllers`), `Test/`, and `ch3_params.m` — the single source of
truth for every knob.

### Why the code order differs from the chapter's

The chapter presents optimization (3) before I/O linearization (4). The code
cannot: the collocation transcription uses `u_ff`, the feedforward produced by
the I/O linearization, as its input. So stage 4 is built first and stage 3
depends on it. `ch3_main` runs them in dependency order.

---

## Design decisions worth knowing

**The control-affine split is exact, not finite-differenced.** The single-
support KKT system has the same left-hand matrix for every input, so
`ch3_control_affine` does one factorization with five right-hand sides and
recovers `ddq = ddq_drift + ddq_in·u` and `λ = lam_drift + lam_in·u` exactly.
This is what makes `L_gL_fy` exact and, in stage 8, lets the friction cone and
the minimum normal force be genuine linear constraints on `u` rather than
approximations.

**The phase variable is linear in q.** `θ = q_t + q₁ + q₂/2 = c·q` is the
absolute stance leg angle. The geometric `atan2` form is *exactly* this
expression for every physical pose (both links are 0.5 m) but wraps at ±π; the
linear form is the same function with the branch cut removed, and its gradient
is a constant row. That constant `ds/dq` is why `L_f²y` has only one curvature
term.

**The phase clamp has a margin.** It engages outside `[−0.5, 1.5]`, not at
`[0,1]`. Clamping hard at the endpoints zeroes `ds/dq` at `s = 0` — the start
of *every* step, where floating point puts `s` on either side of zero — which
silently drops the `(dy_d/ds)ṡ` term from `ẏ`. A polynomial evaluated a whisker
outside `[0,1]` is perfectly well behaved.

**Bézier vs B-spline.** A clamped B-spline of degree M with M+1 control points
*is* the degree-M Bézier curve, and `ch3_test_vc` asserts both values and
derivatives agree. Bézier is the default because its derivatives are
closed-form and it provides an **analytic second derivative**; the B-spline
path central-differences `d²y_d/ds²`, which `L_f²y` depends on directly.

> Building this package surfaced a latent bug in the inherited
> `BSpline_derivative.m`: correct at degree 3 (the only degree the existing
> pipeline calls) but **wrong at degree 5**, drifting up to 0.4 from a finite
> difference of its own curve. Its recursion loop range shrank with degree and
> silently dropped the top basis functions. Fixed, with
> `Test/test_bspline_derivative.m` validating degrees 2–5.

**The collocation solve runs under pure feedforward.** The gait is designed
*on* the zero dynamics surface `Z = {η = 0}`, where any feedback term
multiplies zero. `p.controller` selects what *runs* the gait, not what designs
it.

**The constrained CLF-QP must be run as sampled data, not continuous
feedback.** Its `u(x)` is only *piecewise* smooth — the QP's active set changes
as torque bounds engage and disengage, and `u` kinks at every switch. An
adaptive explicit solver reads each kink as a failed error test and shrinks its
step without bound. This does not merely slow the simulation, it stalls it:
measured, **51 908 RHS evaluations advanced 0.0015 s of a 0.3009 s step** (0.5%),
with `h` collapsed to ~3e-8.

Setting `p.control_dt > 0` solves the QP once per control period and holds it,
so within a period the integrand is the smooth `f + g·u` with `u` constant.
The same step then completes in **1.6 s**. This is also the more faithful model
— the chapter's own framing is a QP solved "well above 1 kHz", which is a
sampled controller, not a continuous-time law. `ch3_compare_controllers` samples
*all three* controllers at 1 kHz, since a continuously-evaluated controller
would otherwise enjoy an advantage no digital implementation of it has.

`ch3_test_simulation` validates the hold against the continuous rollout of the
same controller, and checks the property that actually distinguishes a correct
zero-order hold from one that is merely close: the error must be **first order
in the period**. Measured halving ratios are 1.99 and 1.99.

**Each controller is judged against its own certificate.** The stage-7
min-norm law is built from the CARE `P` and satisfies that rate by
construction — measured, it rides its bound at a ratio of exactly 1.0000. PD
does *not* inherit that rate; its matching certificate is the Lyapunov
equation `AᵀP + PA = −Q` with `A = [0 I; −K_p −K_d]`, which is the form the
chapter writes. Measured against the CARE certificate instead, PD transiently
exceeds the bound by 3.5× before converging further overall. Pairing them the
wrong way is a category error, not a bug.

---

## The trap: small defects do not mean a real trajectory

**Read `ch3_col_verify` before trusting any collocation result.**

Small Hermite–Simpson defects mean the *discrete* equations are satisfied. They
do not mean the nodes approximate a solution of the ODE. On a mesh too coarse
for the dynamics, the optimizer will happily find a **spurious discrete
solution**. Measured here, on a fully converged N = 15 solve:

```
interval-1 defect                      7.18e-07
|x₂(collocation) − x₂(true flow)|      3.68e-02      ← five orders worse
```

The tell was `max|η| = 0.41` at the interior nodes when node 1 satisfied
`η = 0` to 5e-7. Since `u_ff` makes `ÿ = 0` *exactly*, `η` can only drift
through discretization error — so large `η` with a small node-1 residual is a
mesh diagnostic, not a modelling error. The accelerations reached 81 rad/s²,
so `dq` moved by ~2.3 within one `h = 0.072 s` interval; Hermite–Simpson's
truncation error had no chance. The forward simulation disagreed with the
collocation (T: 1.00 s vs 0.72 s), which is how it surfaced.

`ch3_report` runs this check automatically and prints **REJECT** when it
fails. The cure is mesh refinement — `ch3_col_remesh` moves a solution to a
finer mesh, warm-starting from the coarse one.

## Periodic is not stable

The periodicity equality makes the start state a *fixed point* of the
step-to-step map. A fixed point can repel. `ch3_poincare` returns the spectral
radius ρ of that map's Jacobian: **ρ < 1 attracts (walks), ρ ≥ 1 repels
(falls)**, no matter how small the periodicity residual is. It costs 26 step
simulations, which is why it is a post-hoc diagnostic rather than a constraint.

## Table 3.1 limits — measure first

The thesis limits are **ATRIAS** numbers: 63 kg with 50:1 harmonic drives, so
its "|u| ≤ 5 Nm" is *motor* torque, 250 Nm at the joint. RABBIT is ~30 kg and
**direct drive** — `u` here *is* joint torque. Copying the numbers across
produces an infeasible problem and a solver that fails for reasons that look
like bugs.

Every limit in `ch3_report` is printed with its **measured** value whether or
not it is enforced, and `[E]` marks the enforced ones. The workflow is: solve
with all of them off, read the measured ranges, then enable them **one at a
time, warm-starting each phase**, in this order:

> **GRF → friction → torque → impulse**

GRF comes first because friction is `|F_x|/F_z`: enabling it while `F_z` still
crosses zero sets the optimizer chasing a division by zero.

---

## The reference gait

`Results/ch3_reference_gait.mat` (`z_opt`, `p`, `R`, `C`) — a periodic,
mesh-verified, **stable** gait, found with NEC1 disabled:

| quantity | value |
|---|---|
| walking speed | 1.174 m/s |
| step length / duration | 0.353 m / 0.301 s |
| periodicity through Δ | 8.1e-07 |
| mesh verification | 6.8e-04 (tol 1e-03) ✓ |
| max \|η\| over nodes | 2.1e-04 |
| stance-foot drift | 8.7e-07 m |
| **Poincaré ρ** | **0.760 → stable** |

The forward simulation reproduces the collocation exactly — step length 0.353
and duration 0.301 on all six steps — which is the cross-check that matters.
All four Table 3.1 quantities pass *without being enforced*: torque 109 Nm,
impulse 12.8 Ns, friction 0.194, min GRF 166 N.

> **NEC1 was the blocker.** Pinning `L_step/T_step = v_des` while the gait shape
> was still far from periodic over-constrained the problem, and the solve
> stalled on spurious discrete solutions. Dropping it and letting the speed
> float converged to a genuine trajectory within 120 iterations. Speed is then
> recovered by continuation, not imposed from the start.

## The stage-8 payoff, measured

`ch3_compare_controllers` runs that one gait under all three laws at 1 kHz,
with a torque box of 65.5 Nm — deliberately *below* the 109 Nm the gait needs,
so it genuinely binds:

| controller | peak \|u\| | over box | max \|η\| | max δ |
|---|---|---|---|---|
| `iolin_pd` | 109.6 | 44.2 | 1.77e-02 | 0 |
| `clfqp` | 203.3 | 137.8 | 2.77 | 0 |
| `clfqp_con` | **65.5** | **0.0** | 17.4 | 71.5 |

Only the constrained QP holds the box, and it does so *by construction*. Two
things are worth reading carefully:

- **The unconstrained CLF-QP is the worst of the three here.** It asks for the
  least-norm `μ` that certifies the rate at each instant, and rides that bound
  at a ratio of exactly 1.0000. But the RES-CLF inequality is a *floor* on
  convergence, not a target, and on this robot that floor is far too slow:
  `c₃/ε = 0.732` is a time constant of **1.37 s against a step of 0.301 s**, so
  meeting it exactly contracts `V` by only 20% per step — while the impact
  expands `V` by up to **34×**. The hybrid budget is `0.80 × 34.1 = 27.3 > 1`,
  so `η` grows step over step and the demanded torque nearly doubles. The
  surplus convergence that min-norm so efficiently eliminates is exactly what
  was paying for stability across the impact.

  **Sampling is not the cause.** Measured on the transverse dynamics, `dt` from
  1 kHz to 100 Hz gives identical rollouts (peak `|μ| = 3.0`, `V(T)/V(0) =
  0.802` at every rate) and `V` never exceeds its certified envelope. The
  continuous-phase guarantee is not violated anywhere — it simply says nothing
  about Δ, and the robot is a hybrid system. Minimum-norm is not the same as
  well-behaved.
- **δ = 71.5 is the point, not a wart.** It is the controller reporting that it
  could not meet the convergence rate inside the actuator limit. Per Remark 3.2
  the exponential guarantee holds only while δ = 0, so this is the theory's
  boundary being crossed visibly — `max |η|` growing to 17.4 is the price.
  A PD law has no comparable mechanism: it saturates and voids its guarantee
  silently.
