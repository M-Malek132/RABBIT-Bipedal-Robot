# Trajectory_Optimization/

**The main working pipeline.** Finds a dynamically feasible, periodic walking
gait by optimizing B-spline virtual-constraint coefficients (and the initial
state) with `fmincon` SQP, following the Hybrid Zero Dynamics method.

## Entry points

| file | role |
|------|------|
| `rabbit_hzd_trajectory_optimization.m` | **Run this.** Single-gait solve: setup → initial guess (+ optional warm start) → `fmincon` → `inspect_solution` → save → animate. |
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
- **Inequality (`c`, 3):** no ground penetration, mid-step swing-foot clearance
  (NIC3), step-length floor.

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
