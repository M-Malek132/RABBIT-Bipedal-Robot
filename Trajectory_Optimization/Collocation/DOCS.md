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
| `col_poincare.m` | **orbital-stability diagnostic**: `[rho,eig] = col_poincare(col_file, controller)`. Full-order step-to-step spectral radius via `poincare_stability` under the chosen controller (`'pd'` fast, `'clfqp'` slow). `rho < 1` ⇒ walkable. |

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
loop. `col_crosscheck` measures the gap. Measured on two N=20 gaits — the plain
solve (`col_result_2026-07-22_22-40-29.mat`) and the one re-solved with the GRF
and impact-impulse toggles ON (`col_result_2026-07-23_00-45-04.mat`,
`p_over = struct('enforce_grf',true,'enforce_impulse',true)`):

| | plain PD | plain CLF-QP | GRF+imp PD | GRF+imp CLF-QP |
|---|---|---|---|---|
| joint-pos RMS err (rad) | 0.037 | 0.032 | **0.026** | **0.017** |
| closed-loop periodicity defect (∞) | 23.5 | 8.9 | 22.3 | 17.3 |

Two findings:

1. **The GRF/impulse constraints make the gait easier to TRACK.** The smoother
   gait roughly halves the joint RMS error under both controllers and makes the
   closed-loop step length nearly match the collocation (vs a ~25% overshoot
   before). CLF-QP tracks tighter than PD throughout.

2. **They do NOT close the orbit.** The periodicity defect stays ~O(20) — a
   post-impact *velocity* mismatch — and does not improve. An earlier note in
   this file blamed a "near-singular impact (~1e4 Ns)"; that was a stale figure
   from an older gait. The actual impulse here was 16.7 Ns (9.8 after the
   re-solve), and making it gentler did **not** help periodicity. The gap is a
   **closed-loop invariance/stability** issue, not an impact-magnitude one:
   direct collocation finds an *open-loop* periodic orbit, but nothing here makes
   that orbit a *stable fixed point of the step-to-step map under a specific
   controller*.

**To actually get a closed-loop-periodic gait** you need a stability condition,
not more physical constraints: add the step-to-step Poincaré spectral-radius
requirement (`poincare_stability`, the `p.enforce_stability` path the shooting
solver already has) in a warm-started final phase, and/or enforce hybrid
invariance of the controlled dynamics. That is the open next step.

> Caveat: `col_crosscheck` compares ONE step from the collocation's node-1
> state. True closed-loop periodicity is a step-to-step Poincaré fixed point
> (`poincare_stability`), a stronger test than the single-step defect here, so
> the absolute defect values above are indicative, not a stability certificate.

## Orbital stability (col_poincare)

`col_poincare(col_file, controller)` extracts the coeffs + node-1 state and
computes the full-order step-to-step spectral radius `rho` via
`poincare_stability` (which now dispatches on `p.controller`). Measured under PD:

| gait | rho | verdict |
|---|---|---|
| plain N=20 (`..._22-40-29`) | **9.83** | unstable |
| GRF+impulse N=20 (`..._00-45-04`) | **19.83** | unstable |

`rho < 1` is required for a walkable gait; both are ~10–20, i.e. the step map
amplifies perturbations ~10–20× per step. This is the quantitative root cause of
the large `col_crosscheck` periodicity defects, and it confirms the point above:
**direct collocation optimizes feasibility with no stability term, so it lands on
an arbitrary (here, violently unstable) member of the periodic-orbit family.**
Note the GRF/impulse constraints made it *more* unstable (9.8→19.8) — physical
realizability and orbital stability trade off.

### Getting a stable gait — the remaining work

Two routes, neither yet run to completion here:

1. **Shooting `p.enforce_stability` (proven, expensive).** Warm-start the
   shooting solver with the stability constraint (`rho <= rho_max`). Each
   fmincon iteration nests `poincare_stability` (~26 sims) inside finite
   differences over 38 vars (~13 min/iter), and driving `rho` from ~10–20 to <1
   is a large change — realistically hours.

2. **Reduced zero-dynamics restricted Poincaré map in the collocation
   (thesis-faithful, cheap, TODO).** Integrate only the 1-DOF underactuated
   dynamics and constrain its scalar return-map eigenvalue. A first attempt
   validated the manifold `θ̈` at interior nodes but hit a decoupling-matrix
   singularity at the impact boundary (s=1); doing it robustly means projecting
   onto the angular-momentum / unactuated direction rather than inverting the
   decoupling matrix. This is the right long-term path (keeps collocation
   ODE-free) but is unfinished.
