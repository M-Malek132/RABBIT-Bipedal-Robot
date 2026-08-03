# Chapter 5 — Safety-Critical Control via Exponential Control Barrier Functions

Chapter 3 built a controller that asks for **stability**. Chapter 4 asked what
survives when the model is wrong. Chapter 5 adds a second, different demand —
**safety**, in the sense of forward invariance of a set — and the whole chapter
is about what that second demand has to look like when the constraint is not
directly actuated.

Chapter 5 deliberately does **not** use RABBIT. The chapter validates the method
on two purpose-built plants and only then, "in the next Chapter", applies it to
3D walking. So `Chapter5/` carries its own two models and shares nothing with
Chapter 3/4 but the code style.

---

## 1. The one idea

Given a constraint function `h(x)`, the safe set is

```
C = { x : h(x) ≥ 0 }                                             (5.1)
```

and the goal is **forward invariance**: start in `C`, stay in `C`, forever. Not
"track a reference that happens to stay inside", and not "penalize leaving" —
stay in it.

The mechanism is the same one that made the CLF useful in Chapter 3: find a
condition that is **affine in the control**, so a quadratic program can carry
it as a row. Section 5.1 finds one. It just does not work here.

---

## 2. Section 5.1 — reciprocal CBFs, and why they are not enough

The reciprocal candidate (5.6) is `B = 1/h`, blowing up at the boundary, with
condition (5.3):

```
Ḃ(x,u) ≤ γ / B(x)
```

Substituting `B = 1/h`, so `Ḃ = -ḣ/h²` and `γ/B = γh`, this is the much more
readable

```
ḣ ≥ -γ h³                     (implemented in this form; see ch5_barrier)
```

`ḣ = L_f h + L_g h · u`, so **the row constrains `u` only through `L_g h`**.

### The failure, as data rather than as a remark

For a constraint of relative degree `r_b > 1`, `L_g h ≡ 0` — that *is* the
definition of relative degree. The barrier row degenerates to

```
0 · u  ≥  -γ h³ - L_f h
```

which contains no control at all. Two outcomes, and
[`ch5_ctrl_cbf_clf_qp`](../Chapter5/Control/ch5_ctrl_cbf_clf_qp.m) distinguishes
them rather than lumping both into "it failed":

| RHS | what happens | why it matters |
|---|---|---|
| ≤ 0 | row is **vacuously true**; the QP solves and returns the plain CLF-QP control | the quiet failure — the barrier does nothing right up until the state crosses the boundary |
| > 0 | row is **unsatisfiable**; the QP is infeasible | the chapter's "or become infeasible" |

Both Chapter-5 plants put their real constraint at relative degree 6 and 4, so
this is not an edge case — it is the generic case. `ch5_test_qp` asserts
`max|L_g h| == 0` on both plants across sampled states.

There is a **second** limitation, less often stated: `B = 1/h` is a barrier only
on `Int(C) = {h > 0}`. At `h ≤ 0` the reciprocal is negative or infinite and
(5.2)'s sandwich between class-K functions has already failed, so a reciprocal
CBF has **no recovery behaviour** if the state is ever perturbed out of the set.
An ECBF condition is a polynomial in the derivatives of `h` and stays perfectly
well defined at `h < 0`. Both facts are checked in `ch5_test_qp`.

---

## 3. Section 5.2 — Exponential CBFs

### 3.1 Virtual Input-Output Linearization

The fix is to stop differentiating once. Differentiate `h` all the way to `r_b`,
which is *by definition* where `u` first appears, and apply the (5.21) split:

```
h^(r_b) = L_f^r_b h + L_g L_f^(r_b-1) h · u
        = L_f^r_b h + L_g L_f^(r_b-1) h · (u_ff + A⁻¹ μ)
        =: b₀(x) + L_b(x) μ  =:  μ_b                              (VIOL)
```

`μ_b` is **affine in μ** — one scalar row, exactly like the CLF row. The barrier
dynamics are then a chain of `r_b` integrators in controllable canonical form,
eq (5.23), and *nothing about the plant survives into it*:

```
η̇_b = F_b η_b + G_b μ_b,   η_b = [h; ḣ; …; h^(r_b-1)],   h = C_b η_b
```

### 3.2 The step from driving to bounding

This is the whole construction, and it is one character wide. To drive `h` to
zero you would apply

```
μ_b = -K_b η_b   ⟹   η̇_b  =  A_b η_b,   A_b = F_b - G_b K_b
```

giving `h(t) = C_b exp(A_b t) η_b(0) → 0`. To keep `h` **non-negative** you
instead *require*

```
μ_b ≥ -K_b η_b   ⟹   η̇_b  ≥  A_b η_b                            (5.25)
```

giving `h(t) ≥ C_b exp(A_b t) η_b(0) ≥ 0`. The equality that was a control law
becomes an inequality that is a constraint, and the exponential that was a
convergence rate becomes a **lower envelope** on the barrier. Hence the name.

`K_b` comes from ordinary pole placement
([`ch5_ecbf_gain`](../Chapter5/Control/ch5_ecbf_gain.m)), on a linear chain,
whatever the plant was. That is what makes Section 5.2 a *design method* rather
than an existence result.

### 3.3 The conditions, and which ones bite

| condition | where | enforced by |
|---|---|---|
| `A_b` Hurwitz | Thm 5.2 | `ch5_ecbf_gain`, error |
| `A_b` **total negative** (real negative eigenvalues) | Thm 5.2 | `ch5_ecbf_gain`, error |
| `x₀ ∈ C_i` for **every** `i = 0…r_b` | Thm 5.1 | `ch5_ecbf_admissible`, checked at the real `x₀` |
| `p_i ≥ -ẏ_{i-1}(x₀)/y_{i-1}(x₀)` | Cor 5.2 | same |

**Total negative is not pedantry.** Theorem 5.1 is proved by applying
Proposition 5.1 `r_b` times down the chain `C_r_b → … → C_0`, and each link uses
the recursion `y_i = ẏ_{i-1} + p_i y_{i-1}` with `p_i` a *real positive number*.
A complex conjugate pair is perfectly Hurwitz and gives the right characteristic
polynomial, but there is no real `p_i` to build the intermediate sets from, so
the induction has no rungs. `ch5_test_ecbf` checks that complex poles are
*refused*, and that `eig(F_b - G_b K_b) == -poles` exactly — not merely that the
result is stable. (`K_b = fliplr(poly(-p))` and `K_b = poly(-p)` are both
Hurwitz; only one places the poles you asked for.)

**Corollary 5.2 depends on `x₀`,** which the chapter itself flags as a real
limitation of the method. Read the rule for what it is: *if the barrier is
already heading toward the boundary at `x₀`, the poles must be fast enough to
have arrested it in time.* Slow poles are not a conservative choice here, they
are an inadmissible one. `ch5_simulate` checks this before it integrates
anything, because running from an inadmissible `x₀` produces a controller that
enforces its row faithfully at every instant and still lets `h` go negative —
which looks exactly like a bug in the barrier and is not one.

### 3.4 The QP

```
min_{μ,δ}  μᵀμ + p δ²                                             (5.31)
s.t.  V̇(η,μ) + λV(η) ≤ δ                     (CLF, slacked)
      A_c(x) μ ≤ b_c(x)                       (input constraints)
      μ_b ≥ -K_b η_b                          (Exponential CBF, HARD)
      h^(r_b)(x,μ) = μ_b                      (VIOL)
```

**The barrier row is never slacked.** Only the CLF row gets `δ`. When the two
conflict — and on the pendulum they genuinely do — safety wins and tracking
waits. A safety constraint you are willing to violate for a smoother input is
not a safety constraint.

By Remark 5.6, the ECBF row is exactly `y_r_b(x) ≥ 0` for the family (5.28):

```
h^(r_b) + a₁h^(r_b-1) + … + a_r_b h  =  (d/dt+p₁)∘…∘(d/dt+p_r_b) h
```

so `qp.y_rb` lets a run be checked against Theorem 5.1's chain of sets rather
than only against `h ≥ 0`. This identity is asserted in `ch5_test_ecbf` from two
independent constructions.

---

## 4. The two plants

|  | `springmass` (Fig 5.1) | `pendulum` (Fig 5.2) |
|---|---|---|
| dynamics | **linear** | **nonlinear** |
| `nx / nu / ny` | 6 / 1 / 1 | 8 / 2 / 2 |
| relative degree | **6** | **4** |
| constraint | `x₃ ≤ x₃ᵐᵃˣ` | `p₂ʸ ≥ p₂ₘᵢₙ` |
| `L_g L_f^(r-1) y` | `k²/(m₁m₂m₃)` | `D(θ)⁻¹ k / J_m` |

Both have `r·ny == nx`, so `η` is a change of coordinates on the *whole* state
and there are **no zero dynamics** — a convenience of these validation systems
that Chapter 6's walking robot will not share.

In both, **the constrained quantity is the controlled quantity**, so the safety
constraint inherits the full relative degree of the plant. There is no cheap
reformulation in which the barrier is relative degree 1. That is exactly the
situation Section 5.1 cannot address.

### Where the numbers come from

The thesis publishes the ECBF poles (`0.12·[10…15]` and `[5 5.5 6 10]`), the
constraint levels, the link lengths ("unit"), and one response number
(`max(x₃) = 3.27 m` for the CLF-QP baseline). It publishes neither the masses,
the stiffnesses, nor the CLF. Everything not given is chosen and recorded in
[`ch5_system`](../Chapter5/Model/ch5_system.m) next to its reason — for example
`λ = 0.20` on the spring-mass is set *against* the published 3.27 (it reproduces
3.29; 0.18 → 3.36, 0.22 → 3.22).

### The pendulum maneuver is not a small one

`θ₁: -π → +π` with `θ₂` held at 0. Both endpoints are the arm pointing
**straight up**, and they differ by a full turn. The only way through is
`θ₁ = 0`, where an unfolded arm puts the end effector at `p₂ʸ = -2 m`. So the
constraint is not a mild bound near the trajectory — it is *incompatible* with
`θ₂ = 0` over part of the maneuver and **forces the arm to fold against what the
CLF is asking for**. That conflict is what Fig. 5.4's "aggressively move the
links" is describing, and it is why the CLF row needs a slack while the barrier
row must not have one.

---

## 5. Results

Run `ch5_main`. Numbers below are from `Results/ch5_*`.

### Serial spring-mass, relative degree 6

| controller | `x₃ᵐᵃˣ` | `max x₃` | `min h` | verdict | `max\|u\|` | %active |
|---|---|---|---|---|---|---|
| CLF-QP | 3.15 | 3.291 | **-0.141** | VIOLATED | 1.31 N | — |
| ECBF-CLF-QP | 3.15 | 3.008 | +0.142 | SAFE | 1.70 N | 10.5% |
| ECBF-CLF-QP | 3.00 | 2.962 | +0.038 | SAFE | 1.74 N | 17.6% |

The baseline overshoots to 3.291 against the thesis's reported 3.27. And
Fig. 5.3's actual caption claim is reproduced: **varying the constraint while
holding the poles fixed leaves peak force and response speed essentially
unchanged** (1.70 vs 1.74 N).

### Two-link pendulum with elastic actuators, relative degree 4

| controller | `p₂ₘᵢₙ` | `min p₂ʸ` | `min h` | verdict | `p99 \|τ\|` | %active |
|---|---|---|---|---|---|---|
| CLF-QP | -1.0 | -1.983 | **-0.983** | VIOLATED | 27 Nm | — |
| ECBF-CLF-QP | -1.0 | -0.953 | +0.047 | SAFE | 22 Nm | 2.8% |
| ECBF-CLF-QP | -0.5 | -0.464 | +0.036 | SAFE | 23 Nm | 5.4% |

Both ECBF runs are **strictly** safe — no boundary excursion at all, and the QP
is feasible at every one of the 30001 samples. The barrier is active for only a
few percent of the run, which is the interesting part: the poles
`[5 5.5 6 10]` are fast relative to the arm, so the row starts pushing well
before `h` reaches zero and the end effector approaches the limit and turns away
rather than sliding along it. The exponential envelope of Definition 5.1 is
doing exactly what it is for.

The cost is visible in `y₂`: the CLF wants `θ₂ → 0` throughout, and in the ECBF
columns `θ₂` is driven to about `-2.3 rad` mid-maneuver — the arm folding. That
is the hard barrier row winning against the slacked CLF row, and it is why the
maneuver takes longer to finish in the tighter case.

### Three things worth being explicit about

**1. The trajectory is not reproducible; the safety property is.** Once the
barrier row is active this closed loop is a *discrete-time sliding mode* — the
state leaves the boundary between samples and is pushed back at the next one —
and which side of an active-set switch a sample lands on is decided far below
the integration error. Measured: rk4 and ode45 separate by ~0.5 rad within 3 s
on the pendulum, while on the linear spring-mass they agree to 1e-14.

So **no claim here rests on an individual pendulum trajectory**, and the figures
should be read as one realization. What both integrators agree on, and what
`ch5_test_ecbf` asserts for both, is that the safe set is never left. The
pendulum therefore defaults to fixed-step rk4: ode45's adaptive stepping buys
no fidelity in a quantity that is not reproducible anyway.

**2. Sampled-data enforcement is only `O(dt)`.** Forward invariance is a
*continuous-time* claim; the row is imposed at sample instants. Measured on the
pendulum at `p₂ₘᵢₙ = -1.0`:

| `dt` | `min h` |
|---|---|
| 2 ms | -2.4e-3 |
| 1 ms | -1.2e-3 |
| 0.5 ms | **+4.7e-2** |

The default is 0.5 ms — the first rate at which the reproduction actually
*exhibits* invariance rather than nearly does. `ch5_test_ecbf` measures the
scaling rather than assuming it away, and `ch5_report` tags an excursion as
`safe~O(dt)` only when it is within `10·dt`, so a real violation cannot hide
behind the label.

**3. The min-norm torque spike.** A min-norm CLF-QP applies `‖μ‖ = ψ/‖L_gV‖`,
which is unbounded wherever `L_gV` passes near zero. On the pendulum baseline
this happens once, at `t = 0.224 s`, where `‖L_gV‖` dips 140× below its median
and the torque touches 435 Nm for **one sample out of 30001** against a 99th
percentile of 27 Nm. `ch5_report` prints `max|u|` beside `p99|u|` and names the
spike; the figures use a robust shared axis and annotate the true peak rather
than compressing every curve to accommodate it or cropping it silently.

### The Section 5.1 comparison

Same plant, same initial condition, two relative degrees:

| constraint | rel. deg. | controller | result |
|---|---|---|---|
| `ẋ₁ ≤ 0.8` | 1 | CLF-QP | `min h = +0.015` (grazes) |
| `ẋ₁ ≤ 0.8` | 1 | CBF-CLF-QP | `min h = +0.373` — **§5.1 works** |
| `x₃ ≤ 3.00` | 6 | CBF-CLF-QP | `min h = -0.712`, QP **infeasible on 13.5%** of samples |

Both §5.1 failure modes appear in the relative-degree-6 run, in order: the row
is *vacuous* while `-γh³ - L_f h ≤ 0`, and then *unsatisfiable* once the cart
nears the limit with speed. That is the chapter's "high control inputs or become
infeasible", and `ch5_main` reports the measured split rather than asserting one
in advance.

---

## 6. Layout

```
Chapter5/
  ch5_params.m               single source of truth
  ch5_main.m                 both studies + the Section 5.1 comparison
  ch5_gen_pendulum.m         symbolic generator — run once, output committed
  Model/
    ch5_system.m             per-plant block; where the chosen numbers live
    ch5_springmass.m         (5.32)-(5.35)
    ch5_pendulum.m           (5.36)-(5.37), symbolic-safe
    ch5_pend_pv.m            the canonical 12-parameter packing
    ch5_x0.m
    ch5_control_affine.m
    Generated/               exact Lie derivatives, DO NOT EDIT
  Control/
    ch5_io_lin.m             (3.9)/(3.13) at arbitrary relative degree
    ch5_res_clf.m            CLF for arbitrary r  (ch3_res_clf is r = 2 only)
    ch5_clf_eval.m
    ch5_barrier.m            h, η_b, and BOTH barrier forms, kept apart
    ch5_lie_rows.m           exact Lie stack for linear plants
    ch5_pole_gain.m          K_b = fliplr(poly(-p)) — the reversal matters
    ch5_ecbf_gain.m          Thm 5.2's conditions
    ch5_ecbf_admissible.m    Cor 5.2, at the real x₀
    ch5_ctrl_clf_qp.m        the baseline that violates
    ch5_ctrl_cbf_clf_qp.m    Section 5.1, incl. its two failure modes
    ch5_ctrl_ecbf_clf_qp.m   (5.31) — the chapter's result
    ch5_control.m
  Simulation/ ch5_ode_rhs.m ch5_simulate.m
  Analysis/   ch5_report.m ch5_plot_springmass.m ch5_plot_pendulum.m
              ch5_animate_pendulum.m ch5_pend_points.m …
  Test/       ch5_test_all.m + model / barrier / ecbf / qp
```

### On the generated code

The ECBF condition needs `L_f^r_b h` and `L_g L_f^(r_b-1) h` **exactly**. For
`r_b = 4` on a 2-link arm those are fourth derivatives of a trigonometric
expression through an inverted mass matrix, and the obvious alternatives are
both bad: nested finite differences lose ~3 digits per level (by level 4 the
error is `O(1e-3)` on a quantity the controller sits *exactly* on), and
complex-step is exact but does not nest. So it is derived symbolically once and
committed — the same convention the repo already uses for `Dynamics/`.

`ch5_gen_pendulum` calls `ch5_pendulum` with `syms` in place of numbers, so
there is exactly **one** statement of the dynamics in the repository and the
generated derivatives cannot drift from the model the simulation integrates.
`ch5_test_barrier` then checks the result two independent ways: against
finite differences along a tightly-integrated drift trajectory (loose, but
sensitive to a wrong derivative), and against a hand-written Faà di Bruno chain
rule (exact, to 1.4e-15).

---

## 7. Dependencies and solver notes

Optimization Toolbox (`quadprog`), and the Symbolic Math Toolbox **only** to
regenerate `Model/Generated`. No Statistics Toolbox — `ch5_prctile` exists so
that printing a table does not add a third dependency.

**`quadprog` runs `active-set` here, not the `interior-point-convex` that
Chapters 3 and 4 use**, and the difference is not marginal. When the barrier
row is pushing hard, the CLF objective and the ECBF constraint disagree by a
wide margin — the CLF alone would use `‖μ‖ ~ 8` where the barrier row demands
`~93`. Over one pendulum run, interior-point failed on **222 of 30001** samples
with exitflag `-3` ("unbounded") on a problem whose Hessian is positive definite
and therefore *cannot* be unbounded. Re-solving those same QPs:

```
interior-point-convex     0 / 222
active-set              222 / 222,   worst residual 1.4e-11
```

Two supporting fixes came out of the same investigation, both recorded where
they live: constraint rows are normalized by their **coefficients only, never
including the right-hand side** (`ch5_scale_row` — including it flattens an
active barrier row to a gradient of 0.011 and is what interior-point choked
on), and every QP is warm-started at the closed-form min-norm CLF control
(`ch5_min_norm_mu`).

### Running it on a flaky machine

`ch5_main` is **resumable**, which is worth knowing if MATLAB is unreliable
where you run it:

```matlab
F = 'Results/ch5_result_<stamp>.mat';
ch5_main('studies', {'springmass'}, 'only', 1, 'resume', F, 'plot', false)
ch5_main('studies', {'pendulum'},   'only', 2, 'resume', F, 'plot', false)
ch5_main('studies', {},             'resume', F)            % figures + GIF
```

`resume` names the output file and loads it if present; `only` selects
individual runs within a study; anything already in the file is never redone.
Results are checkpointed after every study.

---

## 8. What this hands to Chapter 6

An ECBF is a barrier row that is affine in `μ`, works at any relative degree,
is designed by pole placement, and sits in the same QP as the CLF. Chapter 6
applies exactly that to 3D walking with varied step length/width and to
time-varying stepping stones. The two limitations to carry forward are the ones
Section 5.2.3 closes with and this implementation measures: **pole choice
depends on the initial condition** (`ch5_ecbf_admissible`), and **the whole
construction is model-based** — which is precisely the gap Chapter 4 spent its
length on, and which nothing here addresses.
