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

Alongside stage 3, `HZD/` implements the hybrid zero dynamics itself —
`ch3_zd_point` (the surface `Z` at one phase) and `ch3_zero_dynamics` (the
restricted Poincaré map, `δ_zero`, `V_zero`, `ζ*₂`) — which is what the §6.3.4
stability conditions are stated in. See the NIC/NEC section below.

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

The §6.3.4 gates slot into the same ordering by how much of the trajectory they
can spread a correction over. `phase_mono`, `decoupling` and `swing_clear` are
cheap and usually already satisfied — turn them on early as guards. `liftoff`
and `impact` go **last**, and `impact` needs a march rather than a switch:

> **grf → friction → torque → impulse → phase_mono/decoupling/swing_clear → liftoff → impact**

**Why `impact` is the hardest gate in the pipeline.** Torque, GRF and continuous
friction are all evaluated *along* the step, so the optimizer has every node to
spread a correction over. The impact impulse is a property of **one** state,
`x_N`, filtered through `Δ` — and `x_N` is not free: periodicity ties it to node
1 and the guard equality pins the swing foot to the ground there. The only way
to change the impulse is to change the whole orbit. `ch3_impact_march` steps
`p.limits.mu_s_impact` down a ladder, warm-starting each stage, the way
`ch3_continuation` marches speed and `ch3_posture_march` marches posture.

`mu_s_impact` is a **separate knob from `mu_s`** on purpose. Physically it is the
same floor, but marching `mu_s` itself would silently drag NIC2 along with it.

---

## Section 6.3.4 — the NIC / NEC constraint set

The full constraint set of Westervelt et al. §6.3.4 is implemented as 17 rows in
`ch3_col_constraints`, each individually gated:

| row | constraint | gate | source |
|---|---|---|---|
| 5 | friction cone `\|F_x\| ≤ μ_s F_z` | `friction` | NIC1 |
| 6 | minimum normal force `F_z ≥ F_z_min` | `grf` | NIC2 |
| 9, 10 | swing foot strictly clear; transversal strike | `swing_clear` | NIC3 |
| — | average walking rate (in `ceq`) | `enforce_nec1` | NEC1 |
| 11 | post-impact swing-leg lift-off | `liftoff` | NEC2 |
| 12, 13 | impulse compressive; impulse inside the cone | `impact` | NEC3 |
| 14, 15 | fixed point exists; fixed point is stable | `hzd` | NEC4/NEC5 = (5.79)/(5.80) |
| 16 | `θ` strictly monotonic | `phase_mono` | HH6 |
| 17 | decoupling matrix invertible on `Z` | `decoupling` | HH2 |

**The hypotheses are the purpose of the constraints, not extra rows.** HGW2 is
discharged by rows 6 and 9, HI3 by rows 11–13, HGW6 by rows 3 and 16. HH4/HH5
(hybrid invariance, `Δ(S ∩ Z) ⊂ Z`) gets **no row**: this transcription pins
`y = ẏ = 0` at node 1 and equates `Δ(x_N)` to node 1, so invariance is *implied*.
Adding it again would duplicate rows the periodicity block already spans — the
same rank-deficiency argument that keeps `y = 0` at node 1 only. It is measured
instead (`E.eta_post`, 3.9e-06 on the reference gait), so the implication is
verified rather than assumed.

**Where the book's "NEC" label is loose.** NEC2–NEC5 are filed under nonlinear
*equality* constraints, but every one is written as an inequality — "is
positive", "ζ*₂ > …", "0 < δ²_zero < 1". They are imposed that way here. NEC1 is
the only genuine equality in the list, and the only one in `ceq`.

### The zero dynamics, and why δ²_zero is worth the trouble

`ch3_zero_dynamics` builds the restricted Poincaré map from `α` alone — no
trajectory, no simulation. `ch3_zd_point` reduces the dynamics on `Z` to a scalar
ODE `a(θ)θ̈ + b(θ)θ̇² + c(θ) = 0` by projecting onto the row `w` that annihilates
*both* the actuators and the contact wrench, then an integrating factor
`m = exp(∫b/a)` converts it to the book's `(κ₁, κ₂)` form with `σ = m θ̇`. The step
map is then affine in `ζ = σ²/2`, so its fixed point and eigenvalue are closed
form.

The payoff is a **stability certificate that costs one quadrature instead of 26
step simulations**. On the reference gait the two agree to five digits:

```
delta_zero^2 = 0.75988   (quadrature over alpha)
Poincare rho = 0.75989   (26 forward step simulations)
```

`ch3_test_hzd` checks the reduction against the full 14-state model
(7.8e-10), and checks that `½σ² + V_zero(θ)` is conserved along a real rollout
(8.7e-07) — an invariant that cannot be satisfied by accident if `m` is wrong.

**`enable.hzd` defaults to false, and the cost is why.** Measured on the
reference gait: one constraint evaluation goes from 0.003 s to 0.241 s, so one
fmincon gradient goes from ~1 s to ~77 s. Also note that this transcription
imposes periodicity *directly*, so a converged solve is already at the fixed
point — NEC4/NEC5 are a **check** on the fixed point it found, not the mechanism
that finds one. The book needs them as constraints because its optimization
parametrizes `α` alone and never propagates a state.

### What this found in the reference gait

**NEC3 fails.** The impact impulse is mostly *horizontal* where a walking impact
should be mostly vertical, giving `|I_x|/I_z = 1.35` against `μ_s = 0.4`. The
impact map imposes `J_sw dq⁺ = 0`, "the foot sticks"; at that ratio it would skid
instead. Continuous-phase friction (NIC2) is a comfortable 0.194, so this is
**invisible unless the impulse is checked separately** — exactly why the book
lists NEC3 apart from NIC2.

**And NEC3 is the one constraint in the set that a coarse mesh gets wrong.**

| mesh | `\|I_x\|/I_z` | `I_z` | `ch3_col_verify` |
|---|---|---|---|
| N = 21 | 2.44 | 4.86 Ns | 6.75e-04 ✓ *passes* |
| N = 41 | 1.35 | 8.34 Ns | 3.57e-04 ✓ |

Remeshing alone — with `impact` still **off**, nothing pushing on the impulse —
moved the ratio by 80%. The N=21 gait is not a spurious discrete solution; it
passes verification. But `ch3_col_verify` bounds `max|X_node − X_true|` over the
step, and the impulse is `Λ(q_N)·v_foot(x_N)` evaluated at **one endpoint**, with
`I_z` small enough that a 7e-04 state error swamps it. Every Table 3.1 quantity
is an extremum over the whole step and survives a coarse mesh. This one does not.

The practical consequence: **a ladder calibrated on the coarse number never
becomes active.** A first attempt starting at 2.20 spent two stages with the
constraint inactive — the ratio drifted *up*, 1.104 → 1.261 — and the optimizer
wandered into a region N=41 could not resolve. Refine first, read the ratio, then
ladder from just below it.

NEC2 also only just holds: the trailing foot lifts off at 0.0495 m/s, two orders
below the gait's own 3.88 m/s strike rate.

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
| **δ²_zero** (NEC5) | **0.75996 → stable**, agrees with ρ to 2e-05 |
| ζ*₂ vs `V_max/δ²` (NEC4) | 3.114 vs 0.601 → holds |
| NEC3 impulse `\|I_x\|/I_z` | **2.44 vs μ_s = 0.4 → fails** |

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
  least-norm `μ` that certifies the rate at each instant, which leaves no margin
  — under sampling that error compounds, and the demanded torque nearly doubles.
  Minimum-norm is not the same as well-behaved.
- **δ = 71.5 is the point, not a wart.** It is the controller reporting that it
  could not meet the convergence rate inside the actuator limit. Per Remark 3.2
  the exponential guarantee holds only while δ = 0, so this is the theory's
  boundary being crossed visibly — `max |η|` growing to 17.4 is the price.
  A PD law has no comparable mechanism: it saturates and voids its guarantee
  silently.
