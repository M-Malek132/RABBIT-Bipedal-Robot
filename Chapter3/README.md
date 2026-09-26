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

> **The robot changed on 2026-09-02.** The dynamics were regenerated and the
> model went from 30 kg to 74 kg (`M(1,1)`). The numbers below were re-measured
> on today's model on 2026-09-16, along the reference gait
> `Results/ch3_gait_posture_195.mat`. Anything marked **(30 kg)** predates the
> regeneration and is kept for the mechanism it shows, not for its value. The
> measured account of the whole chapter is
> [`docs/ch3_report.html`](../docs/ch3_report.html) (Persian:
> [`docs/ch3_report_fa.pdf`](../docs/ch3_report_fa.pdf)).

---

## Quick start

```matlab
startup                 % from the repo root; adds Chapter3/ to the path
ch3_test_all            % verify every stage (~6 s)
```

To look at the reference gait, the one verified on today's model and the one
Chapter 4 uses:

```matlab
S = load('Results/ch3_gait_posture_195.mat');
p = ch3_upgrade_params(S.p);
ch3_report(S.z, p);
ch3_plot_gait(S.z, p, 'Results/ch3_gait.png');
```

`out = ch3_main;` solves a gait end to end from the analytic seed and reports on
it, but read this first:

> **No cold start has produced a gait on today's model.** The one cold-start
> campaign on record, `ch3_stage3_from_scratch`, left 18 solves in
> `Chapter3/Results/s3_*.mat` (committed 2026-09-10), and none of them verifies
> as a real trajectory: they crawl at 0.13–0.30 m/s with 1.2–1.5 s steps.
> Everything verified on this model descends from one older gait re-converged
> on the new dynamics — see *The reference gait* below.

> **Long solves and this MATLAB install.** Solves here have been observed to
> die mid-run from crashes inside MATLAB's own add-on registry and worker
> threads — unrelated to this code, but fatal to the process. `ch3_col_solve`
> therefore checkpoints `z` every `p.checkpoint_every` iterations to
> `p.checkpoint_file`. If a run dies, load the checkpoint and carry on; nothing
> is lost but the last few iterations. Running `matlab -nojvm` avoids one of
> the two observed crash paths.

---

## What you call, and what you don't

Chapter 3 is 64 `.m` files, but only **17 are entry points**. The other 47 are
internals reached through them — you should rarely need to call one directly.
If you are looking for "where do I start", it is one of these.

**Solve and verify**

| Call | What it does |
|------|--------------|
| `ch3_main` | solve a gait end to end, then report on it |
| `ch3_test_all` | run all seven test suites (`ch3_test_params`/`_model`/`_vc`/`_control`/`_collocation`/`_simulation`/`_hzd`) |
| `ch3_params` | every knob in the pipeline; the single source of truth |

**Optimization campaigns** — hand-run drivers that march one requirement at a
time, warm-starting each solve from the last. See *Table 3.1 limits — measure
first* below for why marching is not optional, and what order to apply them in.

| Call | Marches |
|------|---------|
| `ch3_continuation` | walking speed `v_des` |
| `ch3_posture_march` | torso-pitch box and hip-height band |
| `ch3_impact_march` | the NEC3 impulse friction cone down to `mu_s` |
| `ch3_lean_tall_march` | combined forward-lean + hip-height campaign |
| `ch3_realizability_march` | finishes Table 3.1: torque, then impulse |
| `ch3_speed_march` | a forward-lean, raised-hip family swept across 0.35–1.20 m/s inside all four Table 3.1 limits |
| `ch3_stage3_from_scratch` | a script, not a function: the staged GRF → friction → torque → impulse solve from the analytic seed |
| `ch3_col_resume` | not a march — runs a bounded chunk of one solve from a checkpoint, one MATLAB process per chunk, for the crash-prone install described above |

The marches warm-start from the seed gaits their headers name —
`ch3_realizability_march` from `ch3_gait_fix`, `ch3_speed_march` from
`ch3_gait_full_constrained` — and neither file is an orbit of today's dynamics
(max `|ceq|` 0.46 and 0.34). Re-converge a seed on the current model before
marching from it.

**Analysis** — the first four take `(z, p)` straight from a solve;
`ch3_animate` takes an unpacked `(x0, alpha, p)`, and `ch3_doc_figures` takes
nothing.

| Call | Produces |
|------|----------|
| `ch3_report(z, p)` | full diagnostic read-out; the first thing to read after a solve (`ch3_main` runs it for you) |
| `ch3_plot_gait(z, p)` | six-panel gait summary figure |
| `ch3_compare_controllers(z, p)` | the stage-8 payoff: same gait under each control law |
| `ch3_hip_accel(z, p)` | hip acceleration from the dynamics, not differenced in time (see below) |
| `ch3_animate(x0, alpha, p)` | walk animation |
| `ch3_doc_figures` | redraws the report's two MATLAB figures into `docs/figures`, at the width `docs/ch3_report_fa.tex` prints them |

Everything else — `Model/`, `VirtualConstraints/`, `Control/`, `HZD/`,
`Simulation/`, the `Optimization/ch3_col_*` transcription, `ch3_seed`,
`ch3_repose`, `ch3_logln`, `ch3_forces`, `ch3_poincare`, `ch3_body_points`,
`ch3_assert_limits`, `ch3_upgrade_params` — is internal. The stage table below says which stage each
one implements.

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

Supporting: `Simulation/` (`ch3_ode_rhs`, `ch3_step`, `ch3_simulate`) runs a
gait forward; `Analysis/` reports on one; `Test/` verifies each stage. Which of
those you call and which are internal is the map above.

Of the analysis tools, `ch3_hip_accel` is the one `ch3_report` does not run for
you: it computes the hip acceleration pointwise from the closed-loop dynamics rather
than by differencing the trajectory in time, because the hybrid motion has a
velocity *jump* at every foot strike and any finite difference across a strike
reports a spike whose height is set by the sample spacing rather than by the
gait. Run it directly when you need that quantity:

```matlab
ch3_hip_accel(out.z_opt, out.p, 4, 'Results/ch3_hip_accel.png');
```

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
measured on the 30 kg model, **51 908 RHS evaluations advanced 0.0015 s of a
0.3009 s step** (0.5%), with `h` collapsed to ~3e-8.

Setting `p.control_dt > 0` solves the QP once per control period and holds it,
so within a period the integrand is the smooth `f + g·u` with `u` constant.
The same step then completes in **1.6 s**. (On today's model
`ch3_test_simulation` runs a constrained-QP step under the hold in 0.3 s.) This is also the more faithful model
— the chapter's own framing is a QP solved "well above 1 kHz", which is a
sampled controller, not a continuous-time law. `ch3_compare_controllers` samples
*all three* controllers at 1 kHz, since a continuously-evaluated controller
would otherwise enjoy an advantage no digital implementation of it has.

`ch3_test_simulation` validates the hold against the continuous rollout of the
same controller, and checks the property that actually distinguishes a correct
zero-order hold from one that is merely close: the error must be **first order
in the period**. Measured halving ratios are 1.97 and 2.00.

**Each controller is judged against its own certificate.** The stage-7
min-norm law is built from the CARE `P` and satisfies that rate by
construction — measured, it rides its bound at a ratio of exactly 1.0000. PD
does *not* inherit that rate; its matching certificate is the Lyapunov
equation `AᵀP + PA = −Q` with `A = [0 I; −K_p −K_d]`, which is the form the
chapter writes. Measured against the CARE certificate instead, PD transiently
exceeds the bound by 2.9× before converging further overall. Pairing them the
wrong way is a category error, not a bug.

---

## The trap: small defects do not mean a real trajectory

**Read `ch3_col_verify` before trusting any collocation result.**

Small Hermite–Simpson defects mean the *discrete* equations are satisfied. They
do not mean the nodes approximate a solution of the ODE. On a mesh too coarse
for the dynamics, the optimizer will happily find a **spurious discrete
solution**. Measured on the 30 kg model, on a fully converged N = 15 solve:

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
its "|u| ≤ 5 Nm" is *motor* torque, 250 Nm at the joint. Copying the numbers
across produces an infeasible problem and a solver that fails for reasons that
look like bugs. RABBIT is not direct drive either: each joint is driven through
a 50:1 harmonic drive and a belt (Chevallereau et al. 2003, Table I and Fig. 4).
In this 74 kg model `u` is the **joint-side** torque, after that reduction; rotor
inertia and gear friction are not modelled (Chapter 4's structured uncertainty
set adds both). The 120 Nm in `ch3_params` is the project's declared limit, not a
Table 3.1 number.

Even the project's own 120 Nm box is out of reach on this model: no verified
gait sits inside it. Two warm-started torque ladders verified down to 212 Nm and
to 195 Nm, then lost verification at their next rungs, 159 Nm and 180 Nm; the
reference gait was solved against its own 195 Nm box. A later ladder from a
1.2 m/s gait reaches a floor of 162.5 Nm, and at that box an impulse ladder
lands a verified gait with a 15.00 Ns impulse (`ch3_impulse_march`; swing apex
0.32 m, stability not measured) — see `docs/ch3_report_fa.tex`.

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

### A gate a saved result never mentioned is **off**

Results carry the `p` that produced them, and `ch3_upgrade_params` fills in
fields added since. `p.limits.enable` is the one nested struct that does **not**
take its missing fields from `ch3_params`: an absent gate is filled in as
`false`. Absence is evidence about the solve — the row was not in the
transcription when that gait was written — not an opinion the result forgot to
record.

This is not hypothetical. `b64160e` added the six NIC/NEC gates defaulting to
`false`, so the merge read older files correctly. `e7e101b` then flipped every
default to `true`, and from that commit the loader switched six constraints on
underneath **eight** stored gaits. `Results/ch3_gait_forward_lean_tall.mat` is
the sharp case: on the 30 kg model it was solved on, it verified as a real
trajectory at `1.30e-05` and missed NEC3 by `0.92` (`|I_x|/I_z = 0.463` against `μ_s = 0.4`), and it is the documented
warm-start seed for `ch3_lean_tall_march`. Nothing was written and nothing
re-solved to cause that — only a default in another file moved.

Two consequences worth keeping:

* **Enabling a gate is always explicit.** `ch3_impact_march` already worked this
  way (it turns `enable.impact` on for its own copy of `p` so the caller's stays
  off); the marches that inherit a seed's gates now say which ones they hold.
* **`ch3_col_verify` cannot catch this.** It asks whether the nodes lie on a real
  trajectory. Being a real trajectory and being feasible for the current
  constraint set are different questions — `ch3_col_check_limits` asks the
  second, and `ch3_assert_limits` refuses to write a deliverable gait that fails
  it.

---

## Section 6.3.4 — the NIC / NEC constraint set

The full constraint set of Westervelt et al. §6.3.4 is implemented in
`ch3_col_constraints`, whose inequality vector has 18 rows (the others hold
Table 3.1 and gait-style limits). Each §6.3.4 row is individually gated:

| row | constraint | gate | source |
|---|---|---|---|
| 5 | friction cone `\|F_x\| ≤ μ_s F_z` | `friction` | NIC2 |
| 6 | minimum normal force `F_z ≥ F_z_min` | `grf` | NIC1 |
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
instead (`E.eta_post`, 7.5e-07 on the reference gait), so the implication is
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
step simulations**. On the reference gait the two agree to 3e-05:

```
delta_zero^2 = 0.74619   (quadrature over alpha, 321 points)
Poincare rho = 0.74621   (26 forward step simulations)
```

`ch3_test_hzd` first checks that the reference gait is still a periodic orbit of
the dynamics on disk (7.6e-07) — a gait file does not record which dynamics it
was solved on, and a stale one would otherwise read as a broken reduction. It
then checks the reduction against the full 14-state model (5.0e-10), and that
`½σ² + V_zero(θ)` is conserved along a real rollout (2.3e-06) — an invariant that
cannot be satisfied by accident if `m` is wrong.

**`enable.hzd` is on by default, and it doubles the cost of an evaluation.**
Measured on the reference gait (N = 61): one constraint evaluation takes 0.017 s
with the gate off and 0.034 s with it on, so one central-difference gradient
over its 879 variables goes from ~30 s to ~59 s. Also note that this
transcription imposes periodicity *directly*, so a converged solve is already at
the fixed point — NEC4/NEC5 are a **check** on the fixed point it found, not the
mechanism that finds one. The book needs them as constraints because its
optimization parametrizes `α` alone and never propagates a state.

### What it measures on the reference gait

On `posture_195` every enforced row holds, and the gait presses on several at
once. Three inequality rows are active to within the solve's 7e-06 tolerance:
peak torque (195.0 Nm against its own 195 Nm box), interior swing clearance
(1.0 mm against 1.0 mm) and the NEC3 impulse cone (`|I_x|/I_z = 0.400` against
`μ_s = 0.4`). The torso pitch sits on the bottom of its box (+1.7°), and the
normal force comes within 0.05 N of its 50 N floor. NIC2 friction is 0.246;
NEC2 lift-off is 0.079 m/s against a 3.18 m/s strike; NEC4/NEC5 hold with
`ζ*₂ = 2.580` against `V_max/δ² = 0.136`. Two limits whose gates were off during
its solve are exceeded: the impulse norm, 15.06 Ns against 15, and the peak swing
height, 0.212 m against 0.15.

The impulse cone being active is the gate doing its job. Continuous-phase
friction says nothing about the single impulse filtered through `Δ` — on the
30 kg model the old reference gait passed NIC2 at 0.194 while its impact ratio
stood at `|I_x|/I_z = 1.35` — which is exactly why the book lists NEC3 apart
from NIC2. That gait also taught the lesson below.

**(30 kg) NEC3 is the one constraint in the set that a coarse mesh gets wrong.**

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

---

## The reference gait

`Results/ch3_gait_posture_195.mat` (`z`, `p`) — a forward-leaning gait that is
periodic, mesh-verified and **stable on today's 74 kg model**:

| quantity | value |
|---|---|
| walking speed | 1.563 m/s |
| step length / duration | 0.427 m / 0.273 s |
| basis / mesh | Bézier degree 5 / N = 61 |
| periodicity through Δ | 4.8e-09 |
| mesh verification | 1.6e-04 (tol 1e-03) ✓ |
| max \|η\| over nodes | 2.7e-05 |
| stance-foot drift | 1.0e-07 m |
| **Poincaré ρ** | **0.746 → stable** |
| **δ²_zero** (NEC5) | **0.746 → stable**, agrees with ρ to 3e-05 |
| ζ*₂ vs `V_max/δ²` (NEC4) | 2.580 vs 0.136 → holds |
| NEC3 impulse `\|I_x\|/I_z` | **0.400 vs μ_s = 0.4 → active** |
| peak torque / impulse | 195.0 Nm (its own box) / 15.06 Ns (gate off) |
| torso pitch / hip height | +1.7 … +4.5° / 0.916 … 0.943 m |

The forward simulation reproduces the collocation — step length 0.427 and
duration 0.273 on all six steps under PD — which is the cross-check that matters.

It was **not** found from scratch. An older 61-node gait, solved before the
regeneration, was re-converged on today's dynamics
(`Results/ch3_gait_fix_reconverged.mat`: verified at 1.3e-05, 1.314 m/s,
287 Nm peak), warm-started down torque ladders to 195 Nm, and then moved through
two posture rungs, torso pitch −5.7…−0.9° → −1.1…+1.7° → +1.7…+4.5°. Every
intermediate is in `Results/` (`gf2_*`, `gf3_*`, `posture_test_result*`).

Like the old reference gait, it was solved with NEC1 off: it carries
`v_des = 0.35` and walks at 1.563 m/s.

> **(30 kg) NEC1 was the blocker.** Pinning `L_step/T_step = v_des` while the gait
> shape was still far from periodic over-constrained the problem, and the solve
> stalled on spurious discrete solutions. Dropping it and letting the speed
> float converged to a genuine trajectory within 120 iterations. Speed is then
> recovered by continuation, not imposed from the start.

## The stage-8 payoff, measured

The experiment `ch3_compare_controllers` runs, re-run step by step on the
reference gait: every law sampled at 1 kHz with the gait's own CLF (CARE,
ε = 0.5), four steps each. Only `clfqp_con` is told about the torque box — 60% of
the gait's 195 Nm peak feedforward, as `ch3_compare_controllers` sets it, and
also 80% and 100%:

| controller | box | steps | max \|η\|, step by step | peak \|u\| |
|---|---|---|---|---|
| `iolin_pd` | not told | 4 of 4 | 0.107, 0.215, 0.214, 0.213 | 195.2 |
| `clfqp` | not told | falls in step 3 | 0.228, 2.44, fall | 283.5 (step 2) |
| `clfqp_con` | 117 (60%) | falls in step 2 | 3.31, fall | 117.0 |
| `clfqp_con` | 156 (80%) | falls in step 2 | 1.60, fall | 156.0 |
| `clfqp_con` | 195 (100%) | falls in step 4 | 0.227, 2.41, 143, fall | 195.0 |

A step counts as a fall when it is shorter than `p.step_len_min` (0.15 m), its
max `|η|` exceeds 1e3, or the state is not finite; steps were capped at 1 s of
robot time (`p.T_max = 0.5`). **The Chapter 3 simulator has no fall detection
of its own** — a guard crossing ends a step even when the robot is on the
ground — so a fallen run keeps integrating to the `2·T_max` cap. With the default
3 s cap, `ch3_compare_controllers` was still inside the unconstrained QP's run
after 40 minutes of wall clock; lower `p.T_max` before running it on a gait that
may fall.

- **The unconstrained CLF-QP is the worst of the three here.** It asks for the
  least-norm `μ` that certifies the rate at each instant, and rides that bound
  at a ratio of exactly 1.0000. But the RES-CLF inequality is a *floor* on
  convergence, not a target, and at ε = 0.5 the floor is slow: `c₃/ε = 0.732` is
  a time constant of **1.37 s against a 0.273 s step**, so meeting it exactly
  contracts `V` by only 18% per step, while every impact kicks `η`. PD, which is
  not tied to that rate, settles at max `|η|` ≈ 0.21 on the same gait. The
  surplus convergence that min-norm so efficiently eliminates is exactly what was
  paying for stability across the impact — which is why Chapter 4 tightens ε
  from 0.5 to 0.20 before it measures anything.

  **(30 kg) Sampling is not the cause.** Measured on the transverse dynamics,
  `dt` from 1 kHz to 100 Hz gave identical rollouts (peak `|μ| = 3.0`,
  `V(T)/V(0) = 0.802` at every rate) and `V` never exceeded its certified
  envelope. The continuous-phase guarantee is not violated anywhere — it simply
  says nothing about Δ, and the robot is a hybrid system. Minimum-norm is not
  the same as well-behaved.
- **The constrained QP does what stage 8 promises, and no more.** The torque
  never leaves the box, by construction, and δ reports the conflict visibly:
  δ = 199 in step 1 at the 60% box, and at 100% it tracks the unconstrained law
  until the box binds, then needs δ = 4.8e5 and 94 QP fallbacks in step 3. Per
  Remark 3.2 the exponential guarantee holds only while δ = 0, so this is the
  theory's boundary being crossed visibly; a saturating PD law voids its
  guarantee silently. But its core is the same min-norm law, and on this 74 kg
  gait no box it was given kept it walking.
