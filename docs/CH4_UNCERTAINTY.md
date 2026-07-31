# Chapter 4 — Robust CLF-QP and L₁ Adaptive Control

Chapter 3 built a controller that assumed it knew the robot. Chapter 4 removes
that assumption and asks what can be done about it. Everything else — the gait,
the virtual constraints, the CLF — is inherited unchanged.

---

## 1. The one idea

Chapter 3 had **one** model. Chapter 4 has **two**, and the whole chapter lives
in the gap between them:

| symbol | name | who uses it |
|---|---|---|
| `f, g` | **true** model | the *simulation* integrates this |
| `f̃, g̃` | **nominal** model | every *controller* is built from this |

This split is enforced structurally, not by convention. `ch4_control` calls
`ch4_io_lin(x, alpha, p, [])` — the empty argument pins it to the nominal
model — and there is no code path from any controller to the perturbed
dynamics. Only [`ch4_ode_rhs`](../Chapter4/Simulation/ch4_ode_rhs.m) touches
both, and its whole body is four lines that make the asymmetry obvious.

Set `p.uncertainty` to zero and Chapter 4 collapses onto Chapter 3 exactly.
That is asserted in `ch4_test_model`, not merely claimed.

---

## 2. What uncertainty does to the Chapter-3 controller

Apply the *nominal* pre-control `u = ũ_ff + (L_g̃L_f̃y)⁻¹ μ` to the *true*
plant, and the clean double integrator `ÿ = μ` becomes **eq (4.3)**:

```
ÿ = μ + Δ₁ + Δ₂μ
```

with (eq 4.4)

```
Δ₁ = L_f²y − L_gL_f y (L_g̃L_f̃ y)⁻¹ L_f̃²y
Δ₂ = L_gL_f y (L_g̃L_f̃ y)⁻¹ − I
```

Two failure modes, and they are different:

- **Δ₁ ≠ 0** — the closed loop no longer has an equilibrium. Tracking error
  cannot reach zero, however hard the controller pushes.
- **Δ₂ ≠ 0** — the controller's *authority* is wrong. At ‖Δ₂‖ ≥ 1 the model can
  cancel or reverse the intended effect, and the feedback can destabilize.

[`ch4_uncertainty`](../Chapter4/Control/ch4_uncertainty.m) computes both.

### The structure that makes this chapter tractable

For a **uniform mass/inertia scale `s`** — the chapter's own perturbation — the
stance KKT system

```
[ sM  −Jᵀ ] [ q̈  ]   [ −sV − sG | B ]
[ J    0  ] [ λ  ] = [ −J̇q̇     | 0 ]
```

splits exactly: substituting `λ = s·w` divides the first row through by `s` and
leaves the second untouched. Therefore

```
q̈_drift  is INDEPENDENT of s          →   L_f²y = L_f̃²y
q̈_in     scales as 1/s                →   L_gL_f y = L_g̃L_f̃ y / s
```

and so, with **no approximation at all**,

```
Δ₂ = (1/s − 1) · I                    exactly isotropic
Δ₁ = −(1/s − 1) · L_f̃²y               exactly proportional to the output drift
```

Verified to ~1e-13 in `ch4_test_model`. Three consequences that shape
everything downstream:

1. **The `'scalar'` Δ₂ model is exact here, not an approximation.** The thesis
   reduces the min–max to two linear inequalities by assuming `Δ₂ = d₂I`
   (Remark 4.5). For uniform scaling that assumption is *true*.
2. **‖Δ₂‖ = |1/s − 1| < 1 requires `s > 0.5`.** Below half the nominal mass no
   worst-case design can help — the bound itself contains a model that cancels
   the control. This is a hard limit of §4.1, not a tuning problem.
3. **Δ₁ tracks the output drift**, which grows away from the orbit. A constant
   `Δ₁max` valid over a *neighborhood* is several times the value valid *on*
   the orbit (measured: ~226 on the gait, ~914 with states jittered by 0.05).
   Since the commanded ‖μ‖ scales as `Δ₁max/(1 − Δ₂max)`, that gap **is** the
   "unnecessarily aggressive" limitation §4.1.4 closes on — with a number on it.

---

## 3. Robust CLF-QP (§4.1)

`V̇` under uncertainty (eq 4.8), using `G = [0; I]`:

```
V̇ = L_fV + L_gV·Δ₁ + L_gV(I + Δ₂)μ
```

The max over the Δ₁ ball is a **constant** — it does not involve μ:

```
max‖Δ₁‖≤D₁  L_gV·Δ₁ = D₁‖L_gV‖
```

so fold it into ψ (Chapter 3's RES residual) and define the **robustified
residual**

```
a := ψ + D₁‖L_gV‖,        ψ = L_fV + (c₃/ε)V
```

The Δ₂ max depends on what the bound means:

| `p.rclf.delta2_model` | worst case | constraint type |
|---|---|---|
| `'scalar'` (default) | `D₂·|L_gV·μ|` | **two linear inequalities**, one per sign |
| `'matrix'` | `D₂‖L_gV‖‖μ‖₂` | second-order cone |

### Both models give the same unconstrained control

The least-norm feasible point lies along `−L_gVᵀ`; substituting `μ = −t·L_gVᵀ/‖L_gV‖`
gives the *same* inequality either way, so

```
μ* = − a·L_gVᵀ / ((1 − D₂)‖L_gV‖²)      when a > 0
μ* = 0                                   when a ≤ 0
```

Setting `D₁ = D₂ = 0` recovers Chapter 3's `μ* = −(ψ/‖L_gV‖²)L_gVᵀ` **exactly**
— the robust controller is a strict generalization, verified to 1e-15.

The `1 − D₂ > 0` requirement is the precise form of Remark 4.3's "within this
region". [`ch4_delta_bounds`](../Chapter4/Control/ch4_delta_bounds.m) reports it;
`ch4_ctrl_rclf_qp` returns `qp.feasible = false` rather than dividing by a
non-positive number.

### Constrained form (4.13) and Remark 4.4

Same structure as Chapter 3 stage 8 — decision vector `[u; δ]`, penalized slack,
torque/friction/GRF rows — **but the added constraints are evaluated on the
nominal model**, so:

- the **torque box IS** invariant to model uncertainty (u is computed from the
  nominal model and applied verbatim), so bounding it is exact;
- the **friction cone and GRF minimum are NOT**. They constrain a *prediction*.
  `ch4_forces` reports `grf_pred_error` — how far the controller's own force
  prediction was from the true one — so the size of that gap is visible rather
  than assumed away. Chapter 8 is where it gets fixed.

`qp.robust_constraints` records which rows are actually robust.

---

## 4. L₁ adaptive control (§4.2)

The robust controller pays its worst-case price *unconditionally*. L₁ instead
**estimates** the uncertainty and cancels it, so a perfect model costs nothing.

Applied input splits as `μ = μ₁ + μ₂`:

- **μ₁** follows a reference model — here Chapter 3's CLF-QP itself, so the
  reference model is *nonlinear with no closed form*. That is the contribution
  of §4.2.2 and why the inherited guarantee is the RES rate, not pole placement.
- **μ₂ = −C(s)θ̂** cancels the estimated uncertainty (eq 4.23).

### The controller has state

Unlike everything in Chapter 3, this is a dynamical system, not static feedback.
[`ch4_l1_state`](../Chapter4/Control/ch4_l1_state.m) defines the 5·n_y = 20
entries once:

| block | size | role |
|---|---|---|
| `eta_hat` | 2n_y | state predictor (4.19) |
| `alpha_hat` | n_y | estimate of the ‖η‖-proportional part of θ |
| `beta_hat` | n_y | estimate of the constant part of θ |
| `mu2` | n_y | filter output — the adaptive control applied |

**Two QP solves per call, and they are not interchangeable.** `μ₁` uses the real
η; `μ̂₁` uses the predictor state (eq 4.21/4.37). Substituting one for the other
would fold the reference model's own tracking behavior into the prediction error
and the adaptation would chase it.

### Why the pieces are what they are

- **Predictor** exists because `η̃ = η̂ − η` is *measurable* while `α̃, β̃` are not.
  Subtracting predictor from plant gives eq (4.24) exactly.
- **Adaptation laws (4.30)** are not free choices — they are precisely the `y`
  that cancels the cross term in the composite Lyapunov function (4.27).
  `Gᵀ P_ε η̃` is recovered as `L_gV(η̃)ᵀ/2` from `ch3_clf_eval`, so there is one
  definition of `P_ε` in the repo, not two.
- **Projection operator** ([`ch4_proj`](../Chapter4/Control/ch4_proj.m)) confines
  the estimates to a ball while preserving inequality (4.29), which is what makes
  the error bound (4.35) finite. Without it the estimate can drift and the result
  is void, not just untidy.
- **Low-pass filter** is what separates *how fast we estimate* from *how fast we
  act*. Γ = 1e4 makes θ̂ fast and ragged; feeding that to the joints would put
  high-frequency content into the ground reaction force and lift the foot.

### At footstrike

η jumps. `eta_hat` is re-seeded to η⁺ (gated by `p.l1.reset_predictor`) so the
adaptation does not read the impact as a phantom uncertainty. **`alpha_hat` and
`beta_hat` are always carried across** — they describe a property of the robot,
and footstrike does not change the robot. Discarding them every step would
restart estimation at ~3 Hz and nothing would ever accumulate.

### Zero uncertainty ⇒ L₁ *is* the CLF-QP, exactly

With a perfect model the prediction error stays at zero, so the adaptation never
moves off its initial condition and `μ₂ ≡ 0`. `ch4_test_l1` asserts this to
machine precision. This is the sharpest statement of L₁'s advantage over the
robust controller, which pays regardless.

### Torque saturation (§4.2.3) binds on μ₁ only

The thesis says so explicitly, and it has a visible consequence: the realized
torque is `ũ_ff + (L_g̃L_f̃y)⁻¹(μ₁ + μ₂)` and only the μ₁ part was inside the box.
`l1.u_box_excess` reports how far past, per call, rather than leaving it to be
discovered from a plot.

---

## 5. Running it

```bash
matlab -batch "ch4_main"
```

Pipeline:

| stage | function | section |
|---|---|---|
| 0 | `ch4_load_gait` | — |
| 1 | `ch4_delta_bounds` — **measure before designing** | (4.4), (4.10) |
| 2 | `ch4_compare_controllers(..., 'robust')` | §4.1.4 |
| 3 | `ch4_compare_controllers(..., 'l1')` | §4.2.4 |
| 4 | `ch4_plot_uncertainty` | Figs 4.2–4.10 |

Stage 1 is not optional and comes first: a robust controller run outside its own
bound is not a robust controller, it is an aggressive one.

```matlab
ch4_main('presets', {'l1'}, 'n_steps', 5)
ch4_main('rclf.delta2_model', 'matrix')     % nested fields take dotted names
ch4_main('l1.omega_c', 100, 'l1.Gamma', 1e5)
```

### Reading the comparison table

The claim is **not** "the robust/adaptive controller is better on average" — on
Case I it need not be, and Remark 4.7 says as much. The claim is that its
convergence behavior is **unchanged across perturbations** while the baseline's
degrades. So read **down** each controller across scales, not across controllers
within a scale. `Vend/Vmx` makes that explicit.

---

## 5a. What this implementation actually measures

Reference gait: `ch3_gait_upright.mat`, T = 0.368 s, v = 1.08 m/s. Three steps
per run, ε = 0.35, control at 1 kHz, `Δ₁max = 250`, `Δ₂max = 0.45`.

### §4.1.4 — the robust CLF-QP reproduces Remark 4.6

| controller | scale 1.0 | scale 1.5 | scale 0.7 |
|---|---|---|---|
| A `clfqp` (min-norm) | 3 steps, ‖η‖→0.484 | **1 step**, →3.36 | **0 steps** |
| B `clfqp_con` (box) | 3 steps, →1.21 | 3 steps, →11.3 | 3 steps, →14.3 |
| **C `rclfqp_con`** | 3 steps, **→0.152** | 3 steps, **→0.115** | 3 steps, **→0.227** |

Step times under C: **0.368, 0.370, 0.378** (scale 1); **0.368, 0.370, 0.366**
(1.5); **0.368, 0.371, 0.369** (0.7) — against a nominal 0.368. The robust
controller's convergence behavior is *flat across the perturbations* while both
baselines degrade by one to two orders of magnitude or fail outright. That is
Remark 4.6, reproduced.

Note also `more robustness costs more mu` in `ch4_test_rclf`: ‖μ‖ = 0.57, 78.8,
209, 470 as the bounds are scaled up from zero. The robust controller pays the
worst-case price *even at zero model error* — §4.1.4's stated limitation, with
a number on it.

### §4.2.4 — L₁ improves on the baseline but does not preserve the rate

| controller | scale 1.0 | scale 0.7 | scale 1.5 |
|---|---|---|---|
| A `clfqp` | 3 steps, →0.484 | 0 steps | 1 step, →3.36 |
| B `l1` | 3 steps, →0.304 | 2 steps, →2.38 | 3 steps, →5.75 |
| C `l1_con` | 1 step, →0.466 | 1 step, →0.435 | 1 step, →1.53 |

L₁ beats the baseline at every perturbation (3 steps vs 1 at scale 1.5; 2 vs 0
at 0.7; lower final error at 1.0), and the estimator demonstrably works —
sampled mid-run at scale 1.5, ‖θ̂‖ vs ‖θ_true‖ reads 22.2/22.6, 41.2/38.4,
87.0/85.9, 90.1/92.3. But it does not hold the convergence rate across
perturbations the way the robust controller does.

**This is a property of the setup, not a defect in the adaptation**, and the
reason is worth stating: L₁'s reference model *is* the Chapter-3 min-norm
CLF-QP. It promises to make the perturbed system behave like that controller —
so when that controller is itself impact-marginal on this gait, L₁ faithfully
reproduces marginal behavior. Diagnostics confirm the failure is at footstrike,
not in the continuous phase: tracking holds at ‖η‖ ≈ 0.15–0.55 through a step,
then the impact throws it to 2.4 and there is not enough convergence rate left
to recover before the next one.

Two things follow. First, ε matters more here than in Chapter 3 — hence
`p.eps = 0.35` rather than the inherited 0.5 (see the note in `ch4_params`).
Second, if you want L₁ to hold its rate under large perturbation, strengthen the
*reference model*: lower ε further, or use the constrained CLF-QP as μ₁.
Tuning Γ or the projection bounds does not help — measured, `alpha_max` from 200
down to 50 moved final ‖η‖ from 9.750 to 9.795.

### Known rough edges

- `rclfqp_con` reports `QPfail` counts (101 of ~2000 samples at scale 1.0,
  box 80). The slack variable should always make the QP feasible, so these are
  `quadprog` exit-flag failures on a badly scaled instance, not genuine
  infeasibility — the fallback saturates the feedforward and the run still
  converges. Surfaced through `F.qp_infeasible` rather than hidden.
- Required friction reaches μ ≈ 0.73 at scale 1.5, above the `p.limits.mu_s`
  of 0.4. The foot would slip on a real surface. Friction is not enforced by
  default; enabling it constrains the *nominal* prediction only (Remark 4.4).

---

## 6. Gotchas

**`ch4_forces` is far more expensive per point than the simulation.** It
re-solves the controller at every sample — two QPs and several KKT
factorizations apiece. A 3-step run at 1 kHz emits ~11 000 solver points;
analysing all of them costs minutes per run and tells you nothing extra, since
the control was *held* at 1 kHz. The `max_samples` argument (default 2000)
decimates; `F.t` and `F.x` report the grid actually used. **Use `F.t`, not
`sim.t`, to index anything F returns.**

**Chapter 4 defaults to sampled-data control** (`p.control_dt = 1e-3`), unlike
Chapter 3. Two independent reasons: the constrained QPs are only piecewise
smooth and stall an adaptive solver (the Chapter-3 lesson), and with Γ ~ 1e4 the
adaptation is far stiffer than the robot, so a continuous run makes ode45 resolve
the *estimator* at every step of the *plant*.

**The L₁ state is advanced by RK4, not Euler** ([`ch4_l1_advance`](../Chapter4/Control/ch4_l1_advance.m)).
At Γ = 1e4 the coupled estimator loop runs at a frequency comparable to the
sample rate itself; Euler there can add energy and diverge, and the divergence
would look like "L₁ is unstable" rather than like an integration artifact.

**Ground reaction forces must come from the true model.** `ch3_forces` computes
λ from the controller's own `aux`, which is correct when there is one model.
Here it would report the *predicted* force — off by 50% at mass scale 1.5.
`ch4_forces` computes λ from the true model and returns the nominal prediction
separately as `lambda_nom`.

**A run that leaves its uncertainty set is not evidence for the robust method.**
`ch4_report` flags `bound_violated`; results from such a run should not be cited
as validating the guarantee, because the controller was outside its own
hypothesis.

---

## 7. File map

```
Chapter4/
  ch4_params.m              all knobs; nested fields take dotted overrides
  ch4_load_gait.m           Chapter-3 gait + Chapter-4 params, consistently
  ch4_main.m                the whole study
  Model/
    ch4_control_affine.m    true OR nominal dynamics; unc=[] defers to ch3
    ch4_impact.m            reset map on the true model
  Control/
    ch4_io_lin.m            I/O linearization of either model
    ch4_uncertainty.m       Δ₁, Δ₂ of eq (4.4)
    ch4_delta_bounds.m      measure them → the bounds of (4.10)
    ch4_ctrl_rclf_qp.m      §4.1, eq (4.12) and (4.13)
    ch4_ctrl_l1.m           §4.2, the control law
    ch4_l1_state.m          the 20-entry controller state, defined once
    ch4_l1_deriv.m          its derivative — the four coupled pieces
    ch4_l1_advance.m        RK4 over one control period
    ch4_proj.m              projection operator of (4.26)
    ch4_control.m           single dispatch point, nominal model always
  Simulation/
    ch4_ode_rhs.m           true plant + nominal controller, augmented state
    ch4_step.m              one hybrid step, with controller state
    ch4_simulate.m          chains steps; per-step random load (Fig 4.11a)
    ch4_is_stateful.m       does p.controller carry state?
  Analysis/
    ch4_forces.m            torques, TRUE forces, estimator signals
    ch4_report.m            one controller vs one perturbed model
    ch4_compare_controllers.m   the §4.1.4 / §4.2.4 sweeps
    ch4_plot_uncertainty.m  Figs 4.2, 4.3, 4.4, 4.6, 4.8, 4.9, 4.10
  Test/
    ch4_test_model.m        true-vs-nominal split and the Δ terms
    ch4_test_rclf.m         the robust guarantee, sampled over the ball
    ch4_test_l1.m           projection, error dynamics, filter, behavior
    ch4_test_all.m
```
