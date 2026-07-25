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
| `col_plot_result.m` | **six-panel report** of one gait: `col_plot_result(col_file, outpng)`. Virtual constraints vs the B-spline, node torques against the torque box, contact force, zero-dynamics phase portrait, swing clearance, stick figure. Prints exitflag / `max|ceq|` / `max c` / speed / peak torque / GRF range. |
| `col_plot_grf.m` | **GRF + friction-cone audit** across gaits: `col_plot_grf(files, labels, outpng)`. Plots `Fz(s)`, `Fx(s)`, the ratio `|Fx|/Fz` against `mu_max` on a log axis, and the friction cone in the `(Fz,Fx)` plane (nodes inside the wedge `o`, slipping nodes `x`). Prints which Table-3.1 toggles each gait was solved with and whether it satisfies them. |

## Usage

```matlab
col_verify_seed();                 % sanity-check the transcription on a seed
rabbit_hzd_collocation(20);        % full solve, 20 nodes, newest PD gait as seed

% re-solve with the Table-3.1 physical toggles, warm-started from a converged gait
rabbit_hzd_collocation(20, 'col_result_2026-07-22_22-40-29.mat', [], ...
    struct('enforce_torque',true,'enforce_friction',true));

col_plot_result();                 % six-panel report of the newest gait
col_plot_grf({'col_result_A.mat','col_result_B.mat'}, {'plain','constrained'}, ...
             'col_grf_compare.png');       % GRF + friction-cone audit
col_poincare('col_result_B.mat','pd');     % orbital stability (rho < 1 = walkable)
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

## Physical realizability — friction cone + torque box

The plain solve leaves the Table-3.1 toggles OFF, and the resulting gait is
**not realizable**: `col_plot_grf` on `col_result_2026-07-22_22-40-29.mat`
reports a worst-case friction ratio `|Fx|/Fz = 6.64` (against `mu_max = 0.4`)
and a peak torque of 174 N·m (against `u_max = 120`). Re-solving with both
constraints ON, warm-started from that gait —

```matlab
rabbit_hzd_collocation(20, 'col_result_2026-07-22_22-40-29.mat', [], ...
    struct('enforce_torque',true,'enforce_friction',true));
```

— converges in 47 iterations (`col_result_2026-07-25_12-01-35.mat`):

| | plain | **friction+torque ON** |
|---|---|---|
| exitflag / `max|ceq|` | 1 / 2.6e-9 | 1 / **7.6e-10** |
| cost `∫u²/L` | 6198 | **2531** |
| max `|Fx/Fz|` (`mu_max=0.4`) | 6.64 ✗ | **0.389** ✓ |
| nodes outside friction cone | 2 / 20 | **0 / 20** |
| peak `|torque|` (`u_max=120`) | 174 N·m ✗ | **108 N·m** ✓ |
| `Fz` min / max | 15 / 882 N | 78 / 632 N |
| `T` / `L` / `v` | 0.999 / 0.355 / 0.355 | 1.044 / 0.371 / 0.355 |

Figures: `Results/col_grf_compare.png` (GRF audit, both gaits) and
`Results/col_report_friction_torque.png` (six-panel report).

Three things worth noting:

1. **The violation is a boundary effect, not a distributed one.** The plain
   gait's friction ratio sits at 0.01–0.1 for most of the step and only blows
   up over the last two nodes, where `Fz` collapses to 15 N while `Fx` is still
   −353 N. Only 2 of 20 nodes are outside the cone.

2. **The cost went DOWN 59% when constraints were added** (6198 → 2531). That
   is backwards for a constrained optimum and is diagnostic: the plain solve
   was sitting in a much worse local minimum, and the friction/torque
   constraints knocked it into a better basin. The plain gait was not merely
   unrealizable, it was also burning ~2.4× the control effort. Peak normal
   force also drops from 2.8 to 2.0 × body weight (314 N), so touchdown is
   gentler even though no impulse constraint was active.

3. **Swing clearance is now pinned to its floor.** Peak clearance falls
   0.364 → 0.118 m and mid-swing sits at exactly 0.0200 m =
   `p.swing_clearance_min`, i.e. that constraint is binding too. Numerically
   fine, but 2 cm of toe clearance is thin for anything but flat ground —
   raise `swing_clearance_min` if you want margin.

Note this run had `enforce_grf` and `enforce_impulse` still OFF; the full
Table-3.1 set has not been solved in one shot.

## Orbital stability (col_poincare)

`col_poincare(col_file, controller)` extracts the coeffs + node-1 state and
computes the full-order step-to-step spectral radius `rho` via
`poincare_stability` (which now dispatches on `p.controller`). Measured under PD:

| gait | rho | verdict |
|---|---|---|
| plain N=20 (`..._22-40-29`) | **9.83** | unstable |
| GRF+impulse N=20 (`..._00-45-04`) | **19.83** | unstable |
| **friction+torque N=20** (`..._12-01-35`) | **4.12** | unstable, but 2.4× better |

`rho < 1` is required for a walkable gait; all three are ≫ 1, i.e. the step map
amplifies perturbations 4–20× per step. This is the quantitative root cause of
the large `col_crosscheck` periodicity defects, and it confirms the point above:
**direct collocation optimizes feasibility with no stability term, so it lands on
an arbitrary (here, unstable) member of the periodic-orbit family.**

**Which physical constraints you add matters — they do not all push the same
way.** The GRF/impulse pair made the gait *more* unstable (9.83 → 19.83), but
the friction/torque pair made it substantially *less* so (9.83 → 4.12). Both
are Table-3.1 realizability constraints, so there is no general trade-off
between physical realizability and orbital stability; an earlier version of
this file claimed one on the strength of the GRF/impulse run alone, and the
friction/torque run is a counterexample from the same family.

The mode structure changes too, not just the magnitude. Top-5 `|eig|`:

```
plain            : [9.829, 9.829, 0.681, 0.254, 0.099]
friction+torque  : [4.121, 1.481, 0.444, 0.444, 0.062]
```

The plain gait has two eigenvalues of *equal* magnitude 9.83 — consistent with
a complex-conjugate pair, i.e. an oscillatory divergence. The friction+torque
gait instead has two distinct real unstable modes, and the second (1.48) is
only marginally outside the unit circle, so a stability constraint would
mostly have to attack the single dominant 4.12 mode. Combined with the shorter
distance to close (4.12 vs 9.83 or 19.83), this gait is the best available
warm start for the stability work below.

> Caveat: `col_poincare` finite-differences the step map under a *specific*
> controller; all three numbers above are PD. `col_crosscheck` found CLF-QP
> tracks tighter, so its `rho` may differ. The CLF-QP run is slow (a QP per ODE
> step × ~26 sims) and has not been done.

### Getting a stable gait — the remaining work

Two routes, neither yet run to completion here:

1. **Shooting `p.enforce_stability` (proven, expensive).** Warm-start the
   shooting solver with the stability constraint (`rho <= rho_max`). Each
   fmincon iteration nests `poincare_stability` (~26 sims) inside finite
   differences over 38 vars (~13 min/iter). Start from the friction+torque
   gait (`rho = 4.12`) rather than the plain one — it is the closest known
   orbit to the unit circle, so the constraint has less distance to close.

2. **Reduced zero-dynamics restricted Poincaré map in the collocation
   (thesis-faithful, cheap, TODO).** Integrate only the 1-DOF underactuated
   dynamics and constrain its scalar return-map eigenvalue. A first attempt
   validated the manifold `θ̈` at interior nodes but hit a decoupling-matrix
   singularity at the impact boundary (s=1); doing it robustly means projecting
   onto the angular-momentum / unactuated direction rather than inverting the
   decoupling matrix. This is the right long-term path (keeps collocation
   ODE-free) but is unfinished.
