# Trajectory_Optimization/

**The main working pipeline.** Finds a dynamically feasible, periodic walking
gait by optimizing B-spline virtual-constraint coefficients (and the initial
state) with `fmincon` SQP, following the Hybrid Zero Dynamics method.

Two solvers for the same gait problem: this directory's **single-shooting**
solver (integrates a step with `ode45` inside each fmincon call; the initial
state is a decision variable), and the **direct-collocation** solver in
[`Collocation/`](Collocation/DOCS.md) (dynamics as Hermite–Simpson defect
constraints, no ODE in the loop; the Bézier/B-spline coefficients are the gait
parameters and the initial state is *not* a variable — the thesis method).

## Entry points

| file | role |
|------|------|
| `rabbit_hzd_trajectory_optimization.m` | **Run this (shooting).** Single-gait solve: setup → initial guess (+ optional warm start) → `fmincon` → `inspect_solution` → save → animate. |
| `Collocation/rabbit_hzd_collocation.m` | **Run this (collocation).** Direct-collocation solve seeded from a PD gait; saves `Results/col_result_*.mat`. See [Collocation/DOCS.md](Collocation/DOCS.md). |
| `build_gait_library.m` | **Speed continuation.** Marches the NEC1 target speed `p.v_des` outward from a seed gait in small steps, warm-starting each solve from the previous converged gait; saves a `gait_library` struct. Produces the raw data behind `Results/gait_library.png`. |

## Problem definition (shared by both entry points)

| file | signature | role |
|------|-----------|------|
| `hzd_problem_setup.m` | `[p, lb, ub, options] = hzd_problem_setup()` | **single source of truth** for parameters `p`, decision-variable bounds, and `fmincon` options (SQP, central finite differences, `ScaleProblem` off). |
| `hzd_cost.m` | `J = hzd_cost(z, p)` | objective: integral of squared torque **per unit step length** (Westervelt eq. 6.43). Returns a large finite value for a failed/degenerate step. |
| `hzd_constraints.m` | `[c, ceq] = hzd_constraints(z, p)` | all nonlinear constraints (see below). |
| `unpack_z.m` | `[coeffs, x_start] = unpack_z(z, p)` | splits the decision vector `z = [coeffs(:); x_start]` (size `nu·n_coeffs + 14 = 38`). |
| `make_coeffs.m` | `coeffs = make_coeffs(x, p)` | initial B-spline control points: interpolate from `x`'s actuated joints to a periodicity-motivated touchdown pose (`relabel_state`). Clamped B-spline ⇒ `spline(s=0)` equals `x`'s joints exactly. |

### Constraints (`hzd_constraints`)

- **Equality (`ceq`, 18):** one-step periodicity (13), stance foot on ground (1),
  stance foot zero velocity (2), `θ_end = θ_plus` (1), **NEC1** walking rate
  `L_step / T_step = v_des` (1).
- **Inequality (`c`, 8):** no ground penetration, mid-step swing-foot clearance
  (NIC3), step-length floor, orbital stability (gated), and the four physical
  constraints below (each gated). Disabled ones are held at `-1` so `c` keeps
  length 8.

### Physical realizability — Table 3.1 (teaching note)

Four constraints (thesis Table 3.1) make a gait physically realizable, each with
its **own toggle** so they can be enabled one at a time:

| toggle | constraint | limit (RABBIT, ~30 kg, no gearing) |
|---|---|---|
| `p.enforce_torque`   | `max\|u\| ≤ p.u_max`            | 120 Nm |
| `p.enforce_friction` | `max\|Fx/Fz\| ≤ p.mu_max`       | 0.4 |
| `p.enforce_grf`      | `min Fz ≥ p.grf_min`            | ≥ 20 N (target ~95) |
| `p.enforce_impulse`  | `\|Fe\| ≤ p.impulse_max`        | 30 Ns |

(swing clearance `p.swing_clearance_min`: raise 0.02 → **0.1 m** for the thesis value.)

### Swing-phase controller — `p.controller`

`hzd_problem_setup` sets `p.controller = 'pd'` (fixed-gain virtual-constraint PD).
Set it to `'clfqp'` and the **whole optimizer runs under the CLF quadratic
program** (`Controller/rabbit_clf_qp_controller.m`): every step sim, the cost,
the constraints, and `gait_forces` all use it, because they route through the
dispatcher `Dynamics/hzd_ode_rhs`. Because the QP enforces the torque box (and
optionally `p.clf_enforce_friction`) *inside* the feedback, a gait optimized
this way is torque-feasible by construction. Cost: a QP is solved at every ODE
step, so a run is far slower than PD — **warm-start from a PD gait, don't
cold-start**. Tunables: `p.clf_eps` (rate), `p.clf_R`, `p.clf_relax`,
`p.clf_enforce_torque`, `p.clf_enforce_friction`.

**Values do NOT transfer from the thesis.** Those are ATRIAS (63 kg, 50:1 gears)
— its `|u|≤5 Nm` is *motor* torque = 250 Nm at the joint. RABBIT is ~30 kg with
no gearing (`u` is the joint torque). Measure first (`inspect_solution` prints
all four), then set limits.

**Order matters — a real diagnosis.** On the reference gait the measured values
were: peak `|u|`≈257 Nm, mean `Fz`≈+309 N (≈ body weight, correct sign), but
**min `Fz`≈−193 N** (the foot lifts!), friction ratio ≈255, impulse ≈9.8 Ns
(fine). The huge friction ratio is **not** slip — `|Fx/Fz|` blows up wherever
`Fz→0`, and here `Fz` goes negative. So **GRF is the root** (`Fz>0`, foot in
compression); friction only becomes meaningful once it holds (`gait_forces`
samples the ratio only where `Fz>1 N`). Enabling friction first made `fmincon`
chase a divide-by-zero and diverge; **enable GRF first**, then friction, then
torque/impulse — warm-starting each phase.

> **Caveat.** The reference gait violates almost all of these by large margins
> and is strongly unstable (`rho ≫ 1`), so incremental repair from that seed may
> not converge. A better initial gait, or the reduced 2-D zero-dynamics
> formulation, is the likely path to a periodic + stable + realizable gait.

## B-spline evaluation

| file | signature | role |
|------|-----------|------|
| `bspline_eval.m` | `[b, db] = bspline_eval(c, s, p)` | evaluate one B-spline output and its `s`-derivative at phase `s`. |
| `BSpline.m` | `N = BSpline(n, p, u)` | Cox–de Boor basis functions (clamped knot vector). |
| `BSpline_derivative.m` | `dN = BSpline_derivative(n, p, u)` | derivative of the basis functions. |

## Diagnostics / helpers

| file | signature | role |
|------|-----------|------|
| `inspect_solution.m` | `report = inspect_solution(z, p)` | prints and returns a full breakdown of a solution: θ sweep, step length/duration, **walking speed**, cost, every constraint residual, worst violation, **stance-foot drift**, and the **stability `rho`** (see below). Use it on any `z_opt`. |
| `poincare_stability.m` | `[rho, eigvals, A] = poincare_stability(coeffs, x_start, p)` | orbital stability of a gait: builds the step-to-step **Poincaré map** Jacobian `A` (13×13, central differences on states 2:14) and returns its **spectral radius `rho`**. `rho < 1` ⇒ attracting/**walkable**; `rho ≥ 1` ⇒ periodic but the robot **falls after a few steps**. ~26 step sims per call. |
| `gait_forces.m` | `F = gait_forces(coeffs, x_start, p)` | replays a gait and returns its **physical quantities** (Table 3.1): peak joint torque, min/mean vertical GRF, max friction ratio, impact impulse. Recomputes `tau`/`lambda` per sample with **the same law the sim used** (`hzd_control_and_dynamics` for PD, `rabbit_clf_qp_controller` when `p.controller='clfqp'`) and the impact impulse via `rabbit_impact_map`. Diagnostic; also feeds the physical constraints. |
| `settle_initial_state.m` | `[x_settled, hist] = settle_initial_state(coeffs, x0, p, nIter)` | Poincaré/fixed-point iteration on the step map. **Disabled by default** (`n_settle = 0`) — it only converges for a *stable* orbit (`rho < 1`) and otherwise diverges, which is exactly what it did on the (unstable) hand-built guess. |

### Periodic vs. stable — the NEC5 lesson (teaching note)

A gait that satisfies the periodicity equality is a **fixed point** of the
step-to-step map `P: x_start → reset(impact(simulate(x_start)))` — but a fixed
point can be **unstable**. The reference seed gait has `max|periodicity| ≈
0.007` (a near-perfect fixed point) yet Poincaré spectral radius `rho ≈ 3`:
simulate it open-loop and `|dq|` grows `6 → 6 → 19 → 891 → …` over successive
steps — it falls over.

This is **Westervelt's NEC5**: a periodic orbit must also be *attracting*
(`rho < 1`). It is the gait analogue of the earlier "step in place"
degeneracy — a missing constraint the optimizer happily exploits. Enabling
`p.enforce_stability` adds the inequality `rho ≤ rho_max` (< 1) so the solver
produces gaits that actually walk.

**Cost / workflow.** The stability inequality is **gated** (off by default)
because each evaluation is ~26 step simulations and `fmincon`
finite-differences it, i.e. hundreds–thousands of sims per iteration.
Recommended: solve a periodic gait *without* it (fast), then a **warm-started
final phase** with `p.enforce_stability = true` to drive `rho` below 1.

## Non-code

- `Trajectory.mat` — committed but **not loaded by any code path**; contents
  unverified (`whos('-file','Trajectory.mat')` to inspect).

## Typical loop

```
hzd_problem_setup ─▶ make_coeffs / warm start ─▶ fmincon(hzd_cost, hzd_constraints)
                                                     │            │
                                          simulate_hzd_gait  simulate_hzd_gait
                                                     ▼
                                          inspect_solution ─▶ save .mat ─▶ animate_hzd_result
```
