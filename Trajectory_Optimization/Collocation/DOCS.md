# Trajectory_Optimization/Collocation/

**Direct-collocation** HZD gait design — the thesis method (Hermite–Simpson),
as an alternative to the single-shooting solver one level up
(`rabbit_hzd_trajectory_optimization.m`). The shooting solver integrates a step
with `ode45` inside every fmincon evaluation and promotes the initial state to
a decision variable; this module instead enforces the dynamics as **defect
constraints** (no ODE in the loop) and makes the **B-spline / Bézier
coefficients** the gait parameters — matching the thesis, where the initial
state is *not* an optimization variable.

## Formulation

Decision vector (see `col_pack`/`col_unpack`), for `N` nodes:

```
z = [ X (2·nq × N) ; U (nu × N) ; Lam (2 × N) ; coeffs (nu × n_coeffs) ; T ]
```

node states, node torques, node **stance-foot contact forces** (explicit), the
spline coefficients, and the step duration. Length `20·N + 25`.

**Dynamics** use the plain EOM with the contact force explicit,
`M q̈ + C + G = B u + J_stᵀ λ` (`col_dynamics`), and the trajectory is tied
together by compressed **Hermite–Simpson defects** — so a step needs no
forward integration.

**Constraints** (`col_constraints`):

| kind | what |
|---|---|
| dynamics | Hermite–Simpson defects between nodes |
| contact | stance foot pinned at **every** node (kills the shooting sim's drift) + zero foot velocity at node 1 |
| **virtual constraints** | `y = q(4:7) − yd(s,coeffs) = 0` at every node ⇒ the gait lies on the hybrid **zero-dynamics manifold**, so `coeffs` parametrizes it (thesis) |
| hybrid | terminal swing-foot strike, `θ(q_N)=θ⁺`, and periodicity through `impact→reset` on all coords but the advancing hip-x |
| task | average speed `L_step/T = v_des`, step-length floor, no penetration, mid-step swing clearance |
| physical (gated) | torque box (as bounds on `U`), friction cone / min-GRF on `Lam`, impact impulse — the Table-3.1 toggles inherited from `hzd_problem_setup` |

Balanced by construction: the explicit `Lam` (2N) variables offset the extra
per-node holonomic + virtual equalities, so vars > equalities (well-posed NLP).

## Files

| file | role |
|------|------|
| `rabbit_hzd_collocation.m` | **driver**: `[z,fval,flag,p]=rabbit_hzd_collocation(N, warm_file, max_iter)`. Seeds from a PD gait, solves, saves `Results/col_result_*.mat`. |
| `col_problem_setup.m` | inherits `hzd_problem_setup` and adds `p.N` + the collocation-vector bounds and fmincon options. |
| `col_pack.m` / `col_unpack.m` | decision-vector (de)serialization. |
| `col_dynamics.m` | `[dq; ddq]` from the EOM with explicit contact force. |
| `col_cost.m` | ∫‖u‖²/L_step via the same Hermite–Simpson quadrature (matches `hzd_cost`). |
| `col_constraints.m` | all equalities + inequalities above. |
| `col_seed_from_pd.m` | warm start: simulate one PD step, resample to `N` nodes, recompute node `u`,`λ`. |
| `col_verify_seed.m` | **test**: reports per-block constraint residuals of a PD seed + an FD check; run this before a long solve. |
| `col_crosscheck.m` | **closed-loop cross-check**: `R = col_crosscheck(col_file, controller, do_plot)`. Extracts `coeffs` + node-1 state, simulates one step under `p.controller` (`'pd'`/`'clfqp'`) with a bounded, robust integrator, and reports joint-error RMS/max, step-length/speed/duration collocation-vs-sim, and the closed-loop periodicity defect. |

## Usage

```matlab
col_verify_seed();                 % sanity-check the transcription on a seed
rabbit_hzd_collocation(20);        % full solve, 20 nodes, newest PD gait as seed
```

Verified: from a PD seed (`N=12`) fmincon drives feasibility `55 → 3·10⁻⁹`
(exitflag 1) — a feasible periodic on-manifold gait. Each evaluation is cheap
(no `ode45`), but with `fmincon` + finite-difference Jacobians the wall-clock
still grows with `N`; analytic/sparse Jacobians (or IPOPT) would be the
production step.

## Closed-loop realizability (col_crosscheck)

A collocation gait is periodic on the zero-dynamics manifold *by construction*,
but that does **not** guarantee a given feedback law reproduces it in closed
loop. `col_crosscheck` measures the gap. On the N=20 gait
(`col_result_2026-07-22_22-40-29.mat`, seeded from the natural-speed PD gait):

| | PD | CLF-QP |
|---|---|---|
| joint-pos RMS err | 0.037 rad | **0.032 rad** |
| closed-loop periodicity defect (∞) | 23.5 | **8.9** |

CLF-QP holds the manifold tighter and lands ~2.6× closer to periodic than the
fixed-gain PD law (its exponential-convergence guarantee at work), but **neither
closes the orbit**: this seed-gait family has a near-singular impact (impulse
~1e4) that amplifies any residual tracking error through the reset map. The
practical fix is to re-solve the collocation with the impact-impulse and GRF
toggles ON (`p.enforce_impulse`, `p.enforce_grf`) so the optimizer finds a
gentle-impact gait a controller can actually hold — the toggles are already
inherited from `hzd_problem_setup`.

> Caveat: `col_crosscheck` compares ONE step from the collocation's node-1
> state. True closed-loop periodicity is a step-to-step Poincaré fixed point
> (`poincare_stability`), a stronger test than the single-step defect here.
