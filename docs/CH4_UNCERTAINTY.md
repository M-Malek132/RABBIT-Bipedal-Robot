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
   the orbit (measured on `posture_195`: 232.6 along a nominal rollout, 1140.6
   with the sampled states jittered by σ = 0.05 — 4.9×).
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

**Advisory is not the same as optional.** It is tempting to strip (4.13) down
to the torque box, which reads like the cleaner version of "robust CLF-QP with
torque saturation". Measured on `posture_195`, that is not physical. Once η ≠ 0
the robust μ\* carries a term of fixed magnitude D₁/(1 − D₂) ≈ 574 rad/s² along
−L_gVᵀ, whose direction flips whenever L_gV passes through zero however small η
is, so a 1 kHz loop chatters. Without the normal-force floor, that chatter
pulled the **true** stance foot into the ground for **17 / 38 / 42%** of the
samples in Cases I–III (min Fz −1621 / −3220 / −3543 N). The stance contact is
integrated as a pin, so the simulation walked on regardless. With the rows on
(nominal Fz ≥ 50 N, |Fx| ≤ 0.4 Fz), the true Fz stayed ≥ 23 N in every case
and the robust law still walked all three steps. So the rows stay on, as
`ch3_params` sets them, and the sweep tables print the true `min Fz` next to
the tracking numbers.

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

**The box is sized to the gait, with headroom.** The thesis's 65 Nm is below the
feedforward of every gait in `Results/`, and a μ₁ box that cannot deliver the
feedforward starves the inner QP; the adaptation then chases the tracking
failure instead of the model error and runs away. So `ch4_load_gait` raises
`p.l1.u_max` to **1.25×** the gait's own peak torque (244 Nm on `posture_195`).
The 25% is measured: at mass scale 1.5 a box exactly at the 195 Nm peak ran
away to 7193 Nm and lost the robot in step 3, while 1.25× walked all three
steps at Vend/Vmx 0.065 (1.5× and 2× changed nothing material).

---

## 5. Running it

```bash
matlab -batch "ch4_main"
```

Pipeline:

| stage | function | section |
|---|---|---|
| 0 | `ch4_load_gait` — **refuses a gait that is not an orbit of today's dynamics** | — |
| 1 | `ch4_delta_bounds` — **measure, then adopt** | (4.4), (4.10) |
| 2 | `ch4_compare_controllers(..., 'robust')` | §4.1.4 |
| 3 | `ch4_compare_controllers(..., 'l1')` | §4.2.4 |
| 4 | `ch4_plot_uncertainty` | Figs 4.2–4.10 |
| 5 | `ch4_animate` | — |

Stage 0 stops the run if the gait's own collocation residuals, re-evaluated on
the `M/V/G` currently on the path, exceed `p.verify_tol` (see §6 for why a gait
file cannot be trusted to know which dynamics it was solved on).

Stage 1 is not optional and comes first: a robust controller run outside its own
bound is not a robust controller, it is an aggressive one. Its result is also
what runs. Unless you pass `rclf.delta1_max` / `rclf.delta2_max`, `ch4_main` sets
the bounds to 1.2× the Cases I–III maxima measured along a nominal rollout of
the loaded gait. Case IV (scale 3) is measured and printed but not folded in,
since it would size the robust law for a perturbation the sweeps never apply.

```matlab
ch4_main('presets', {'l1'}, 'n_steps', 5)
ch4_main('rclf.delta2_model', 'matrix')     % nested fields take dotted names
ch4_main('l1.omega_c', 100, 'l1.Gamma', 1e5)
ch4_main('animate', false)                  % skip the GIFs
```

### The animation

[`ch4_animate`](../Chapter4/Analysis/ch4_animate.m) runs several controllers on
**the same perturbed robot** and draws them side by side, synchronized in time.
This is the picture the tables cannot give you: "completed 1 step" and
"completed 3 steps" are two numbers, but one panel frozen on the floor while the
next two keep walking is the actual result.

Synchronizing takes a little work and the details matter. The runs have
different solver grids *and* different durations, so each is interpolated onto
one common clock; a run that has already ended freezes at its last pose, drawn
washed out and labelled `STOPPED at <t>`. Interpolation touches `q` only and
only for drawing — nothing feeds back into a simulation, and an early-ending run
is marked rather than extrapolated. Each panel also reports live `‖η‖`, because
two stick figures look similar for a while before one of them falls, and the
tracking error separates them well before that.

`ch4_main` writes one GIF per mass scale in `anim_scales` (default 1.5 and 0.7),
and hands each one the torque boxes the sweeps used at that scale; runs are
configured by the same `ch4_run_params` as the table rows, so each panel is the
controller its row scored. Standalone, the boxes are whatever `p` carries:

```matlab
[x0, alpha, p] = ch4_load_gait();
p.uncertainty.mass_scale = 1.5;
ch4_animate(x0, alpha, p, {'clfqp','rclfqp_con','l1'}, 4, 'Results/ch4_walk.gif')
```

### Reading the comparison table

The claim is **not** "the robust/adaptive controller is better on average" — on
Case I it need not be, and Remark 4.7 says as much. The claim is that its
convergence behavior is **unchanged across perturbations** while the baseline's
degrades. So read **down** each controller across scales, not across controllers
within a scale. `Vend/Vmx` makes that explicit.

**The robust boxes are measured, not copied.** §4.1.4 sets each case's box
"slightly below the maximum torque that controller A uses" — 60/80/150 Nm on the
thesis's 32 kg robot. On the current 74 kg model those numbers are 12–41% of
what controller A draws, and both constrained laws fell in every case, the
perfect-model one included. The `'robust'` preset applies the *rule* instead:
controller A runs first at each scale, and

```
box = max( 0.8 × peak|u| of A ,  scale × p.gait_u_peak )
```

The floor is the true robot's own feedforward peak (s times the mass needs s
times the torque for the same motion). It binds only in Case I, where A with a
perfect model draws barely more than the feedforward itself (196.4 vs 195.0 Nm):
at 0.8× (157 Nm) `clfqp_con` fell in step 3 and `rclfqp_con` in step 1, while at
the 195 Nm floor both walk all three. Pass `opts.u_box` to fix the boxes by hand.

**`min Fz` is the column that says whether any of it is walking.** It is the
*true* model's normal force. The stance contact is integrated as a pin, so a
controller can demand that the ground pull the foot down and still complete
steps in simulation; a negative entry means those steps are not realizable.

**Read `Vend/Vmx` loosely for `rclfqp_con`.** Its worst-case term chatters at the
sample rate, so V at the final sample lands anywhere in an order-of-magnitude
band. The same Case I run gave 0.021 at a 195.0 Nm box and 0.19 at 195 Nm plus a
floating-point hair. Steps, max‖η‖ and min Fz survive that.

**`steps` counts steps, not guard crossings.** A falling robot reaches the
ground too, and before this was checked a run with max‖η‖ = 290 was tabulated as
"3 steps completed". `ch4_simulate` rejects a crossing — and ends the run there —
when the guard fired the instant it was armed (the swing foot never left the
ground), the hip fell below half height, the torso pitched past 1 rad, or the
step was shorter than `p.step_len_min`. The `reason` field names which.

---

## 5a. What this implementation actually measures

Reference gait: `ch3_gait_posture_195.mat`, the forward-lean 195 Nm gait solved
on the current 74 kg model. T = 0.2733 s, L = 0.427 m, v = 1.56 m/s, and its own
collocation residuals re-evaluated on today's dynamics are 7.6e-7 (defect),
7.6e-9 (periodicity) and 1.2e-6 (‖η⁺‖). Three steps per run, ε = 0.35, control at
1 kHz, `Δ₁max = 279.1` and `Δ₂max = 0.514` (measured along a nominal rollout,
×1.2). Result set `Results/ch4_{robust,l1}_2026-09-12_23-46-32/`, 69 s end to end.

Each cell reads **steps · max‖η‖ · min Fz**. A fall is named by the step
`ch4_simulate` rejected.

### §4.1.4 — the robust CLF-QP reproduces Remark 4.6

| controller | Case I, ×1 (box 195 Nm) | Case II, ×1.5 (box 623 Nm) | Case III, ×0.7 (box 405 Nm) |
|---|---|---|---|
| A `clfqp` (min-norm) | 3 · 0.39 · 39 N | **falls in step 3** · 15.1 · 74 N | **falls in step 3** · 26.7 · −1810 N |
| B `clfqp_con` (box) | 3 · 0.43 · 50 N | **falls in step 3** · 15.1 · 78 N | **falls in step 2** · 3.1 · 12 N |
| **C `rclfqp_con`** | 3 · 1.59 · 50 N | 3 · **0.49** · 63 N | 3 · **2.62** · 23 N |

Both baselines fall under either perturbation. A falls in step 3 when its swing
foot fails to clear at ×1.5, and loses the robot 13 ms into step 3 at ×0.7.
B falls in step 2 at ×0.7 with the hip at 0.17 m. The robust controller walks
every case: step times **0.273, 0.273, 0.273** s (Case I), **0.273, 0.273,
0.274** (II) and **0.274, 0.279, 0.273** (III), against a nominal 0.2733. It
keeps the true normal force positive throughout (≥ 23 N) and max‖η‖ at or below
2.62. Its
convergence behavior is flat across the perturbations while the baselines lose
the robot — that is Remark 4.6.

It is *not* the tightest tracker in Case I (max‖η‖ 1.59 against A's 0.39),
which is Remark 4.7 and the price §4.1.4 closes on. `more robustness costs more
mu` in `ch4_test_rclf` puts a number on it: ‖μ‖ = 0.31, 87.6, 233, 524 as the
bounds are scaled up from zero at one state, paid even at zero model error. On a
1 kHz loop that price shows up as chatter, which is the rapid oscillation in the
yellow V curve of Figure 1 (see §3 for what the contact rows prevent it from
doing).

### §4.2.4 — L₁ keeps the robot stepping, but not on physical contact forces

| controller | ×1 (box 244 Nm) | ×0.7 | ×1.5 |
|---|---|---|---|
| A `clfqp` | 3 · 0.39 · 39 N | **falls in step 3** · 26.7 · −1810 N | **falls in step 3** · 15.1 · 74 N |
| B `l1` | 3 · 0.39 · 49 N | 3 · 21.1 · **−2688 N** | 3 · 10.6 · **−410 N** |
| C `l1_con` | 3 · 0.40 · 50 N | 3 · 21.4 · **−2881 N** | 3 · 9.1 · **−484 N** |

With a perfect model L₁ is the CLF-QP, as §4 says: the same max‖η‖. Its
estimate is nonzero only through the 1 kHz sample-and-hold error, and the lower
Vend/Vmx (0.066 against 0.36) suggests it partly cancels that error. Under either
perturbation it completes three steps at near-nominal timing where its own
reference model falls: 0.271, 0.276, 0.278 s at ×0.7 and 0.278, 0.290, 0.245 s
at ×1.5 for `l1`.

**Those steps are not realizable.** The true normal force is negative for **26% /
15%** of the samples under `l1` (×0.7 / ×1.5) and **25% / 12%** under `l1_con`. μ₂
acts outside any QP, so nothing bounds the contact force it implies, and the
pinned stance model absorbs the pull. Until that is addressed, the §4.2.4 result
here is a tracking result, not a walking one.

Tracking at ×0.7 is also poor (max‖η‖ 21, final ‖η‖ 12–13), and it is not the
box: every box from 1.0× to 2.0× the gait's peak gives max‖η‖ 21–24. It *is*
sensitive to the filter. At ω_c = 50 / 100 / **150** / 300 rad/s, max‖η‖ reads
11.3 / 11.6 / **21.1** / 27.3 at ×0.7 and 7.1 / 10.3 / **10.6** / 9.4 at ×1.5.
150 rad/s is kept because it is the thesis's value and three steps per setting
is thin evidence, but it is the first knob to revisit, and whether a lower ω_c
also removes the negative Fz is not yet measured. The estimator itself follows
the uncertainty's magnitude rather than tracking it pointwise, as §4 expects:
sampled at ×1.5, ‖θ̂‖ against ‖θ_true‖ reads 87/38, 61/71, 58/33, 226/234, 94/70.

L₁'s reference model *is* the Chapter-3 min-norm CLF-QP, and on this gait that
controller falls under both perturbations. Strengthening the reference model
(lower ε, or the constrained CLF-QP as μ₁) remains the structural lever, beyond
tuning the filter.

### Known rough edges

- **Remark 4.4, measured.** The robust law's friction rows hold the *nominal*
  prediction at μ ≤ 0.4, and in Cases I–II the true demand matches (p95 0.40).
  At ×0.7 the lighter true robot needs **p95 1.24, peak 1.40** — the foot would
  slip on a μ = 0.4 floor while the controller believes it will not. That is the
  gap Chapter 8 exists to close.
- `QPfail` counts are small — `rclfqp_con` 2 / 2 / 0 and `l1` 15 / 8 / 7 samples
  of ~2000. The slack keeps the QP feasible in principle, so these are
  `quadprog` exit-flag failures, and the fallback saturates the feedforward.
- `ch4_forces` re-solves each controller at every ODE output point, while the
  simulation holds each torque for 1 ms. For the chattering robust law those
  torques differ pointwise, but not in what the contact columns say: on the
  torque-only runs of §3, recomputing min Fz under the torques actually held
  matched `ch4_forces` to within 53 N on every run.
- Results in `Results/` stamped between 2026-09-02 and 2026-09-12 23:37 predate
  one or more of this section's fixes — a stale gait, copied boxes, falls
  counted as steps, contact rows stripped — and should not be cited.

---

## 6. Gotchas

**A gait file does not record which dynamics it was solved on.** `M.m`, `V.m`
and `G.m` are global, so a gait saved before they were regenerated still loads,
unpacks and simulates — as a reference that is no longer an orbit of the robot.
The 2026-09-02 regeneration took the model from 30 kg to 74 kg, non-uniformly
(torso 10 → 47 kg, femurs 5 → 10, tibias 5 → 3.5), and the then-default
`ch3_gait_upright.mat` went to a collocation defect of 6.6e-2, a periodicity
residual of 1.02 and a verify deviation of 14.1. Chapter 4 ran on it for ten
days: with a *perfect* model the baseline needed 382 Nm and ratcheted ‖η‖ up at
every impact, and the robust sweep collapsed while its table still printed
"3 steps completed". `ch4_load_gait` now re-evaluates the gait on the current
dynamics (`meta.orbit`) and `ch4_main` refuses one that fails. After any
dynamics regeneration, re-solve or at least re-verify the gait first.

**1 kHz sample-and-hold kicks every impact, even with a perfect model.** On
`posture_195`, ‖η⁺‖ after the first step is **0.242 / 0.123 / 0.049** at a
1 / 0.5 / 0.2 ms control period, against 1.7e-4 under continuous control. The
pre-impact error is small (6.8e-3 at 1 ms); the impact map amplifies it about
35×. That is the jump every Case I curve shows at each footstrike — the
discretization, not the model. Lower `p.control_dt` if a figure needs it clean;
runtime grows roughly in proportion.

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
  ch4_load_gait.m           Chapter-3 gait + Chapter-4 params, consistently;
                            re-checks the gait is an orbit of today's dynamics
  ch4_main.m                the whole study
  Model/
    ch4_control_affine.m    true OR nominal dynamics; unc=[] defers to ch3
    ch4_case_dynamics.m     per-scale re-derived M/V/G (0.5, 0.7, 1.5, 3)
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
    ch4_simulate.m          chains steps; rejects falls; per-step random load (Fig 4.11a)
    ch4_is_stateful.m       does p.controller carry state?
  Analysis/
    ch4_forces.m            torques, TRUE forces, estimator signals
    ch4_report.m            one controller vs one perturbed model
    ch4_run_params.m        one comparison run's params: law, perturbation, box, rows
    ch4_compare_controllers.m   the §4.1.4 / §4.2.4 sweeps; measures the robust boxes
    ch4_plot_uncertainty.m  Figs 4.2, 4.3, 4.4, 4.6, 4.8, 4.9, 4.10
    ch4_animate.m           controllers racing on one perturbed robot
  Test/
    ch4_test_model.m        true-vs-nominal split and the Δ terms
    ch4_test_rclf.m         the robust guarantee, sampled over the ball
    ch4_test_l1.m           projection, error dynamics, filter, behavior
    ch4_test_all.m
```
