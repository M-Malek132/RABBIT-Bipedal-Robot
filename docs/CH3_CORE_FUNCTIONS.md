# Chapter 3 — the three core functions

`ch3_control_affine`, `ch3_io_lin` and `ch3_impact` between them contain almost
all of the non-obvious reasoning in the package. Everything else assembles what
they produce. This is a deep read of each: the derivation, the code, the measured
numerical behaviour, and where the comments and the numbers disagree.

Inventory and call graph: [`CH3_FUNCTION_MAP.md`](CH3_FUNCTION_MAP.md).
Stage 3: [`CH3_OPTIMIZATION_FLOW.md`](CH3_OPTIMIZATION_FLOW.md).

All measurements below are on `Results/ch3_reference_gait.mat` unless stated;
the scripts are described in [How the numbers were measured](#how-the-numbers-were-measured).

---

## 1. `ch3_control_affine` — building `f` and `g`

[`Chapter3/Model/ch3_control_affine.m`](../Chapter3/Model/ch3_control_affine.m)

### What it has to produce, and why

Chapter 3 needs the continuous phase as `ẋ = f(x) + g(x)·u`. Not "given `u`,
compute `ẋ`" — that is enough to *simulate* but not to *design*. Three later
stages need the `u`-dependence separated symbolically:

- stage 4 needs `L_gL_f y = Jy·q̈_in`, the sensitivity of `ÿ` to torque;
- stage 7's QP needs `ÿ` affine in `u` to be a QP at all;
- stage 8 needs the ground reaction force affine in `u`, so the friction cone is
  a **linear** inequality rather than an approximation.

### The physics

Unconstrained, the robot obeys `M(q)q̈ + V(q,q̇) + G(q) = B·u`. In single support
the stance foot is pinned, `P_st(q) = const`, which differentiated twice gives
`J_st·q̈ + J̇q̇ = 0`. Holding the foot takes a force, entering as the generalized
force `J_stᵀλ`:

```
M q̈ + V + G = B u + J_stᵀ λ
J q̈         = −J̇q̇
```

Two unknowns, two equations. Stack them into a KKT / saddle-point system:

```
[ M   -Jᵀ ] [ q̈ ]   [ Bu − V − G ]
[ J    0  ] [ λ  ] = [   −J̇q̇     ]
```

`λ` appears as a Lagrange multiplier, so the contact force is not modelled — it
is *solved for*, out of the same linear system as the accelerations.

### The trick

The left-hand matrix contains no `u`, and `u` enters the right-hand side
linearly. So one factorization with five right-hand sides — one drift, four
input directions — recovers the entire affine structure by superposition:

```matlab
rhs = [ [-V_vec - G_vec], B_mat ; ...
        [-Jdotdq],        zeros(nc, nu) ];

sol = A \ rhs;             % (nq+nc) x (1+nu)
```

```
q̈ = q̈_drift + q̈_in·u          λ = λ_drift + λ_in·u
```

Exactly — not linearized, not finite-differenced. This *is* the solution,
rearranged. Verified against [`rabbit_constrained_dynamics`](../Dynamics/rabbit_constrained_dynamics.m)
at **4.5e-13**.

The vector fields then fall out:

```matlab
f = [dq;              ddq_drift];
g = [zeros(nq, nu);   ddq_in   ];
```

The structural zero in the top block of `g` — torque cannot change position
instantaneously — is exactly what makes the outputs relative degree 2 in stage 4.

### Why not finite-difference it

Five solves plus truncation error, and every downstream quantity inherits the
error: `L_gL_f y` → `u_ff` → the collocation cost → every gradient `fmincon`
estimates. Measured:

```
one solve with 5 RHS  9.6 us | five separate solves 33.2 us (3.5x)
```

3.5× faster **and** exact.

### Numerical character

```
gait   : cond(A) min 456.7  max 490.8 | min rank(J_st) = 2
random : cond(A) median 428.2  max 1.073e+03 | rank(J_st) < 2 in 0 of 4000
```

Condition number around 450–500, barely varying, and `J_st` never lost rank.

Two structural notes:

- **`A` is not symmetric** (`|A − Aᵀ| = 5.88`), because the top-right block is
  `−Jᵀ` while the bottom-left is `+J`. Flipping one sign would give the symmetric
  indefinite form that admits `LDLᵀ`. The current arrangement is equally correct
  and yields `λ` directly in the convention you want — **force the ground applies
  to the robot**, `Fz > 0` pushing up (measured `Fz = +39.75 N` at `u = 0`). At
  9×9 the solver choice is irrelevant.
- **The solve is unguarded.** `ch3_io_lin` checks `rcond` and falls back to
  `pinv`; this function checks nothing. The justification is sound — `M` is
  positive definite always and `J_st` for a planar point foot is generically full
  row rank — and 0 rank failures in 4000 random poses supports it. Worth knowing
  it is there.

### Two small things done right

```matlab
M_mat  = M(q);
V_vec  = V([q; dq]);
G_vec  = G(q);
```

The suffixes are not decoration. `M`, `V`, `G` are function names on the path;
`M = M(q)` shadows the function and the next call in scope indexes the matrix
instead.

```matlab
nc = size(J, 1);
```

Contact count read from the Jacobian rather than hardcoded, so the KKT assembly
and the extraction slices follow a different contact model automatically.

`aux` is built only when `nargout > 2`. `ch3_ode_rhs` calls with two outputs and
skips the struct entirely — which matters when `ode45` calls it tens of thousands
of times per step.

### In one sentence

> The KKT matrix does not contain `u`, so one factorization with five right-hand
> sides recovers the exact affine dependence of both the accelerations *and* the
> contact force on torque — and that exactness is what lets stage 4 build an exact
> decoupling matrix and stage 8 write a genuine linear friction cone.

---

## 2. `ch3_io_lin` — from a nonlinear robot to `ÿ = μ`

[`Chapter3/Control/ch3_io_lin.m`](../Chapter3/Control/ch3_io_lin.m)

### What it does

It is the join point: stage 1 and stage 2 come in, the linearized error dynamics
comes out. The whole body of the function is four lines.

```matlab
[~, ~, aux]   = ch3_control_affine(x, p);
[y, ydot, o]  = ch3_outputs(x, alpha, p);

Lf2y  = o.Jy * aux.ddq_drift - o.curv;
LgLfy = o.Jy * aux.ddq_in;
```

### The derivation

The outputs have **relative degree 2**: `u` is a torque, so it cannot appear
until you have differentiated down to `q̈`.

```
y  = H q − y_d(s(q))      → no u
ẏ  = Jy(q) q̇              → no u
ÿ  = J̇y q̇ + Jy q̈          → q̈ carries u
```

Because `ds/dq` is a **constant row** (θ linear in `q` — hypothesis HH6), the
first term collapses to a single curvature term:

```
J̇y   = −(d²y_d/ds²)·ṡ·(ds/dq)
J̇y q̇ = −(d²y_d/ds²)·ṡ²
```

Substituting `q̈ = q̈_drift + q̈_in·u` and grouping:

```
ÿ = L_f²y + (L_gL_f y)·u

  L_f²y    = Jy·q̈_drift − (d²y_d/ds²)ṡ²      (4×1)  what happens at u = 0
  L_gL_f y = Jy·q̈_in                          (4×4)  the decoupling matrix
```

Had θ been nonlinear in `q`, `J̇y` would carry the **Hessian of the phase** — a
7×7 second-derivative tensor to build and differentiate at every timestep. That
is the concrete cash value of "linear in q".

Both quantities are exact: `q̈_in` from the KKT solve, `d²y_d/ds²` analytic from
[`ch3_bezier`](../Chapter3/VirtualConstraints/ch3_bezier.m). No finite differences
anywhere in the chain.

### The inversion

```matlab
u = u_ff + (L_gL_f y)⁻¹ μ,     u_ff = −(L_gL_f y)⁻¹ L_f²y
```

substituted back gives `ÿ = μ`. Verified independently — finite-differencing `ẏ`
along the true flow `f + g·u`, never reusing the `ÿ` formula:

| check | result |
|---|---|
| `u = u_ff` → `ÿ = 0` | max \|ÿ\| = 2.55e-09 |
| `μ = [1.5, −2, 0.7, 3]` → `ÿ = μ` | max error 2.48e-09 |
| `μ = e_k` moves only channel `k` | identity matrix to 1e-9 |

`u_ff` is the torque that *walks*: hold the four constraints, let the pendulum do
the rest. Peak over the reference step: **106.6 Nm**.

Every controller in stages 5–8 differs only in how it picks `μ`. That claim is
enforced structurally by [`ch3_control`](../Chapter3/Control/ch3_control.m), where
all four branches receive the identical `Lf2y, LgLfy, u_ff`.

### What is *not* linearized

On the constraint manifold single support is 10-dimensional. `η = [y; ẏ]` is 8.

```
10 − 8 = 2      ← the zero dynamics, (θ, θ̇)
```

This is **partial** feedback linearization, and partial by necessity: full-state
feedback linearization needs full actuation, and there are 4 actuators for 5
effective DOF. The division of labour for the whole chapter follows:

| | dim | behaviour | handled by |
|---|---|---|---|
| `η = [y; ẏ]` | 8 | linear, `ÿ = μ` | stages 5–8, feedback |
| zero dynamics | 2 | fully nonlinear | stage 3, choice of `alpha` |

### Numerical character, and a correction to the docstring

The guard:

```matlab
rc = rcond(LgLfy);
if ~isfinite(rc) || rc < 1e-12
    u_ff = -pinv(LgLfy) * Lf2y;
```

Measured:

```
gait      : rcond min 5.61e-03  max 8.34e-03 | cond min 80.2  max 128.0
random 3k : rcond median 8.45e-03  min 2.04e-06 | below 1e-12 in 0 of 3000
random 20k: worst rcond = 8.05e-06 | below 1e-12 (pinv path taken): 0
```

**The `pinv` fallback never fires** — not on the gait, not in 20 000 random
states spanning the full ±π joint box. It is genuinely defensive.

The docstring says the decoupling matrix loses rank *"in practice near knee lock,
or if dyd/ds drives Jy toward a rank-deficient combination."* Measured, those two
halves do not fare equally:

**Knee lock: not supported.** Sweeping the stance knee to full extension leaves
conditioning flat —

```
  q2 = 1.347  rcond = 5.762e-03      q2 = 0.100  rcond = 6.829e-03
  q2 = 0.600  rcond = 5.190e-03      q2 = 0.000  rcond = 7.033e-03
```

A 121×121 grid over **both** knees across `(−π, π)` puts the worst case at

```
min rcond = 2.143e-07 at (q2,q4) = (2.125, 1.866)
```

— deep double-knee **flexion** (≈122° and ≈107°), not extension. And still five
orders above the `1e-12` threshold.

**Profile slope: supported.** Stretching `alpha` about `alpha_0` degrades
conditioning monotonically:

```
  scale   1  rcond = 5.762e-03  cond =   128.0
  scale  10  rcond = 2.874e-04  cond =  2668.3
  scale 100  rcond = 1.014e-04  cond =  8252.8
```

A 100× steeper profile costs about 1.75 orders of magnitude of `rcond`. So
`dy_d/ds` really is the term that drives `Jy` toward rank deficiency; the knee
attribution appears to have been a plausible guess rather than a measurement.
**The docstring has been corrected** to state the profile-slope mechanism, to
record that knee angle barely matters, and to note that the `pinv` fallback has
never been observed to fire.

---

## 3. `ch3_impact` — the same structure, applied to an instant

[`Chapter3/Model/ch3_impact.m`](../Chapter3/Model/ch3_impact.m)

### The same matrix, a different meaning

```
[ M   -J_swᵀ ] [ q̇⁺  ]   [ M q̇⁻ ]
[ J_sw   0   ] [ imp ] = [   0   ]
```

Structurally identical to `ch3_control_affine` — same saddle-point shape, same
`M`, a foot Jacobian and a multiplier. Everything it *means* is different:

| | `ch3_control_affine` | `ch3_impact` |
|---|---|---|
| unknowns | accelerations `q̈` | post-impact **velocities** `q̇⁺` |
| multiplier | contact **force** [N] | contact **impulse** [Ns] |
| Jacobian | `J_st`, the pinned foot | `J_sw`, the foot about to land |
| RHS | forces `Bu − V − G` | momentum `M q̇⁻` |
| balance | force | momentum |

It is a momentum balance, not a force balance. `q` does not appear on the
right-hand side because **configuration is continuous through impact** — the
robot does not teleport, only its velocities jump. `M` and `J_sw` are both
evaluated at the pre-impact configuration.

The bottom row imposes `J_sw q̇⁺ = 0`: the foot strikes and **sticks**. No
rebound, no slip — a rigid plastic impact.

Measured: `cond(A) = 490.8`, essentially the same as the continuous KKT system.
`J_sw q̇⁺ = 4.8e-16`.

### It is violently dissipative

```
KE: 70.064 -> 20.269 J   dissipated 49.795 J (71.1%)
impulse: [-11.851  4.858] Ns, norm 12.808
```

**71% of the kinetic energy vanishes in an instant.** That is the physical
content of "plastic": the ground absorbs whatever it takes to stop the foot dead.
This is why the impact is the dominant disturbance in the hybrid loop, and why
stage 5's `ε` is squeezed — the controller has exactly one step to undo what Δ
does in zero time.

The impulse norm of 12.8 Ns is the Table 3.1 impulse quantity.

### Order is not negotiable

Δ is three operations and the docstring insists on the order. It is right, and
the failure mode is worse than "inaccurate":

```matlab
% --- 1. plastic impact ---
sol     = A \ b;
% --- 2. relabel ---
x_rel = ch3_relabel(x_post, p);
% --- 3. re-plant ---
x_plus(2) = x_plus(2) + foot(2);
```

Swap 1 and 2 and here is what happens:

```
wrong-order impulse [0.000 0.000] vs correct [-11.851 4.858]
||dq+ correct - dq+ wrong|| = 13.8518  (correct norm 5.4861)
correct J_sw*dq+ = 4.82e-16 | wrong, on the true swing foot = 0.6631
```

**The impulse is exactly zero — the impact becomes a no-op.** The reason is
worth understanding rather than memorizing: `ch3_relabel` swaps the leg
coordinate slots, so after relabelling, `J_sw` evaluated on the new coordinates
describes the *old stance* foot. That foot was pinned, so it already satisfies
`J·q̇ = 0`. The constraint the solver is asked to enforce is already satisfied,
the multiplier comes out zero, and nothing happens.

The velocity error is 13.85 against a correct norm of 5.49 — a 250% error — and
the foot that actually landed is left sliding at 0.66 m/s instead of stuck. The
code runs, produces plausible-looking numbers, and is silently wrong. Hence:

```matlab
% ORDER MATTERS.  The impact uses the PRE-relabel swing Jacobian J_sw.
% Relabel first and the impulse is applied to the wrong foot.
```

### The re-plant, and the one coordinate left alone

```matlab
foot   = P_st(x_rel(1:nq));
x_plus = x_rel;
x_plus(2) = x_plus(2) + foot(2);
```

`y` is down-positive, so `foot_z = −y − (joint terms)` and adding `foot(2)` to `y`
drives the new stance foot to exactly `z = 0`. On a collocation solution the
correction is negligible — node N is constrained to the guard, so the residual
measured `1.7e-12 m` — but in forward simulation the event tolerance leaves a real
residual, and without this the contact point would ratchet downward across steps.

`px` is deliberately **not** reset:

```
px delta through Delta = 0.00e+00
```

It is the translation gauge and must accumulate so the robot advances. This is
why [`ch3_col_constraints`](../Chapter3/Optimization/ch3_col_constraints.m) writes
periodicity as `E.x_next(2:end) - X(2:end, 1)` — 13 conditions, not 14. The
missing one is the point of walking.

---

## How the numbers were measured

Each figure above comes from a script run against `Results/ch3_reference_gait.mat`:

| claim | method |
|---|---|
| affine split exact | compare `f + g·u` and `λ_drift + λ_in·u` against `rabbit_constrained_dynamics` for random `u` |
| `cond(A)`, `rank(J_st)` | assemble the KKT matrix at every gait node and at 4000 random poses |
| 5-RHS timing | 2000 repetitions of `A\rhs` versus 2000×5 single-column solves |
| `ÿ = μ` | central difference of `ẏ(x)` along `ẋ = f + g·u`, never reusing the `ÿ` formula |
| `rcond(LgLfy)` | gait nodes, 3000 plausible states, 20 000 states over the full ±π box, a 121×121 knee grid, and an `alpha` stretch sweep |
| impact dissipation | `½q̇ᵀMq̇` before and after the KKT solve, `M` at the pre-impact configuration |
| order matters | relabel first, then run the same momentum balance, and compare `q̇⁺` and the impulse |

Model-level invariants behind all of this — mass-matrix symmetry and positive
definiteness, `M(1,1)` = total mass, `DM = dM/dt`, the passivity identity
`q̇ᵀ(Ṁ − 2C)q̇ = 0` at 1.0e-14, `G = ∂P/∂q`, all four Jacobians against finite
differences, and energy conservation at 1.1e-14 over 0.25 s of free rotation about
the toe — were checked separately and all hold.

Two documentation defects found while verifying, neither affecting behaviour:

- `p.mass` in [`ch3_params.m`](../Chapter3/ch3_params.m) read `32` against an
  actual total mass of **30 kg** (10 torso + 4×5 links, confirmed by `M(1,1)`).
  The field is assigned and never read anywhere in the codebase.
  **Corrected to 30.**
- The `ch3_io_lin` knee-lock attribution, discussed above. **Corrected** — the
  docstring now states the profile-slope mechanism the measurements support,
  drops the knee-lock claim they contradict, and records that the `pinv` path
  has never been observed to fire.

Separately, verifying the smoothness of [`ch3_yd`](../Chapter3/VirtualConstraints/ch3_yd.m)
turned up two real defects in the B-spline branch, **both fixed**: the clamp
returned a `dyd` frozen at the boundary slope rather than zero, contradicting the
constant `yd` it accompanied; and the second-derivative stencil was clipped
against `[0,1]`, degenerating to a one-sided formula and costing five orders of
accuracy at `s = 1` (touchdown). See the smoothness section of
`RABBIT_documentation.tex` and `Visualization/make_yd_smoothness_figure.m`.
