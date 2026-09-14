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

### The exact law is sliding-mode control, so it chatters: the boundary layer

Near the orbit ψ is O(‖η‖²) while D₁‖L_gV‖ is O(‖η‖), so the closed form tends to

```
μ*  →  − M · L_gVᵀ / ‖L_gV‖,        M = D₁ / (1 − D₂)  ≈ 574 rad/s²
```

That is a correction of fixed magnitude along a unit vector that flips whenever
the state crosses L_gV = 0, however small η is: unit-vector sliding-mode
control. Held for T = 1 ms it cannot slide. One sample moves L_gV by
φ_T = ‖2P₂₂‖·M·T = 0.029, where P₂₂ is the ẏẏ block of P_ε, and whatever was left
of L_gV is overshot. Over the first three steps of `rclfqp_con` (measured at the
then-default ε = 0.35; nothing in this subsection's algebra depends on ε), the torque
actually held changed by a median of **253 / 573 / 481 Nm per sample** in
Cases I–III.

`p.rclf.boundary_layer` (κ, default 1) replaces the unit vector by a saturation:

```
D₁‖L_gV‖  →  D₁‖L_gV‖ · min(1, ‖L_gV‖/φ),        φ = κ (1 + D₂) φ_T
```

Inside the layer the robust part of μ* is the linear feedback −(M/φ)L_gVᵀ.
One held sample of it multiplies L_gV by 1 − (1 + d₂)φ_T/φ, which is at worst
1 − 1/κ over the Δ₂ bound. That fixes the units of κ:

| κ | sampled loop along L_gV, for every d₂ in the bound |
|---|---|
| 0 | the exact law of (4.12)/(4.13), chatter included |
| > ½ | converges |
| ≥ 1 | converges without overshoot, i.e. without chatter |

The same sweep, measured at the control samples, over the first three steps.
Baseline A reads max‖η‖ 0.39 in Case I and falls in the other two.

| κ | median held-torque change [Nm/sample], I / II / III | max‖η‖, I / II / III |
|---|---|---|
| 0 | 253 / 573 / 481 | 0.81 / 1.16 / 2.54 |
| 0.33 | 1 / 2 / **514** | 0.11 / 0.24 / 3.02 |
| 0.66 | 1 / 2 / 1 | 0.09 / 0.25 / 0.38 |
| 0.99 | 1 / 2 / 1 | 0.11 / 0.42 / 0.48 |
| 1.32 | 1 / 2 / 1 | 0.12 / 0.61 / 0.61 |
| 3.30 | 1 / 2 / 1 | 0.15 / 1.61 / 1.28 |

At κ = 0.33 the multiplier predicts −1.0 / −0.33 / −1.86 for the three true
input gains (1, 1/1.5, 1/0.7), and only the last case, the one below −1,
chattered. A thicker layer trades tracking for margin; κ = 1 is the thinnest
layer that cannot overshoot anywhere in the bound.

Two consequences:

- **Inside the layer the robust correction does not depend on D₁.** It is
  −L_gVᵀ/(κ(1 + D₂)‖2P₂₂‖T). A larger Δ₁ bound only moves the layer's edge out, so
  near the orbit the "more robustness costs more μ" price of §4.1.4 is capped by
  the sample period.
- **The guarantee becomes ultimate boundedness.** (4.11) still holds exactly
  wherever ‖L_gV‖ ≥ φ. Inside the layer the worst case may exceed it by at most
  D₁‖L_gV‖(1 − ‖L_gV‖/φ) ≤ D₁φ/4, reported as `qp.bl_gap`, so
  V̇ ≤ −(c₃/ε)V + D₁φ/4. The thesis's "tracking errors converging to zero" needs
  continuous control. Any implementation that holds its control for a sample
  gets a bounded error instead.

`ch4_test_rclf` checks 9–12 assert each of these: the relaxed guarantee over the
sampled ball; the vanishing correction at the orbit, where the exact law keeps
M; no chatter in a sampled loop that makes the exact law chatter; and the
D₁-independence. `control_dt = 0` also gives φ = 0, since continuous control has
no sample to overshoot.

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
torque saturation". Measured on `posture_195`, that is not physical. With the
exact law (κ = 0) the chatter described above, without the normal-force floor,
pulled the **true** stance foot into the ground for **17 / 38 / 42%** of the
samples in Cases I–III (min Fz −1621 / −3220 / −3543 N). The stance contact is
integrated as a pin, so the simulation walked on regardless. With the rows on
(nominal Fz ≥ 50 N, |Fx| ≤ 0.4 Fz), the true Fz stayed ≥ 23 N in every case.
The boundary layer removes the chatter but not the need for the rows: at the
default layer and ε = 0.20, 25 steps at ×0.7 with the rows off still put 1.8% of
the samples below zero (min Fz −481 N), against ≥ 23 N with them on. So the rows stay on,
as `ch3_params` sets them, and the sweep tables print the true `min Fz` next to
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
| `eta_hat` | 2n_y | state predictor — (4.19), or its plant-input form (below) |
| `alpha_hat` | n_y | estimate of the ‖η‖-proportional part of θ |
| `beta_hat` | n_y | estimate of the constant part of θ |
| `mu2` | n_y | filter output — the adaptive control applied |

### Three fixes to §4.2's formulation

Run as written, §4.2's controller does not reproduce its own claim here. Over 25
steps at ε = 0.20 it tracks no better than its reference model under either
perturbation, and at ×1.5 it falls where the reference model does not. It also
needs the ground to pull the foot down on 11–14% of the samples. Three changes
follow, each an option in `p.l1`. The thesis form stays one setting away:
[`ch4_l1_opts`](../Chapter4/Control/ch4_l1_opts.m) resolves the options, and a
parameter struct saved before they existed resolves to the thesis form.

**1. The predictor is driven by what the plant received** (`p.l1.predictor`).
The thesis predictor (4.19) runs its own copy of the reference model:
`η̂̇ = Fη̂ + Gμ̂₁ + G(μ₂ + θ̂)`, with μ̂₁ a second CLF-QP solved on η̂.
Subtracting the plant gives (4.24), which still contains μ̂₁(η̂) − μ₁(η).

The step to (4.28) needs that difference to decrease P_ε along η̃. The two QPs
only guarantee a decrease along η and along η̂ separately, and the min-norm law
is nonlinear. So the mismatch forces η̃ exactly as a model error would, and the
adaptation absorbs it: θ̂ overshot θ 3–5× at ε = 0.35.

The `'plant'` predictor is the standard L₁ state predictor, with μ = μ₁ + μ₂ as
applied:

```
η̂̇ = Fη + G(μ + θ̂) − a(η̂ − η)
```

Its error obeys `η̃̇ = −aη̃ + Gθ̃` identically, with no reference model in it.
Three things follow:

- **Adaptation signal:** the P in the adaptation law becomes I, so the law
  adapts on the velocity prediction error.
- **Tuning:** the β channel is the loop s² + as + Γ, so a = 2√Γ gives critical
  damping.
- **Cost:** the second QP per call goes away.

**2. The sampled advance sees the plant's inputs** (same option). The thesis
advance froze η across each 1 ms period while letting μ₂ move with the filter.
The plant saw the opposite: η moving, μ₂ held.

Under `'plant'` the advance holds the applied μ and reads η at both ends of the
period. It stays causal, because η is sampled at t_{k+1} before the next control
is computed. Freezing η would inject a·ÿ·τ into the velocity channel.
`ch4_test_l1` check 9 runs a tracking double integrator: reading both ends
estimates θ to 1e-15 relative error, and freezing η leaves 2.5%.

After each period the estimates are also returned to their projection balls.
Projection bounds them only in continuous time; with a stiff loop an RK4 step
can leave the ball, and θ̂ once reached 1e7 before this was added.

**3. The constrained law bounds the torque it applies**
(`p.l1.constrain_applied`). §4.2.3 boxes μ₁ alone and adds no contact rows, so
μ₂ reaches the joints outside every constraint.

With the option on, μ₂ enters the QP as a known offset: the QP's feedforward
becomes `ũ_ff + (L_g̃L_f̃y)⁻¹μ₂`. The cost and the CLF row still act on μ₁,
while the box and the friction and normal-force rows bound the total torque.
The rows are written on the nominal model (Remark 4.4), but they are no longer
bypassed.

Because the box now bounds the total torque, the L₁ preset scales it with the
robot, `u_max · max(1, s)`. A 1.5× robot needs 1.5× the torque for the same
motion, and a 244 Nm box is below the 293 Nm its feedforward alone needs.

§5a has the measured effect.

### Why the pieces are what they are

- **Predictor** exists because `η̃ = η̂ − η` is *measurable* while `α̃, β̃` are not.
  Subtracting predictor from plant gives eq (4.24) for the thesis form and
  `η̃̇ = −aη̃ + Gθ̃` for the plant form, both exactly.
- **Adaptation laws (4.30)** are not free choices — they are precisely the `y`
  that cancels the cross term in the composite Lyapunov function (4.27), with
  the P that certifies the predictor's error. For the thesis form that is P_ε:
  `Gᵀ P_ε η̃` is recovered as `L_gV(η̃)ᵀ/2` from `ch3_clf_eval`, so there is one
  definition of `P_ε` in the repo, not two. For the plant form it is I.
- **Projection operator** ([`ch4_proj`](../Chapter4/Control/ch4_proj.m)) confines
  the estimates to a ball while preserving inequality (4.29), which is what makes
  the error bound (4.35) finite. Without it the estimate can drift and the result
  is void, not just untidy.
- **Low-pass filter** is what separates *how fast we estimate* from *how fast we
  act*. Γ = 1e5 makes θ̂ fast and ragged; feeding that to the joints would put
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

### Torque saturation (§4.2.3), and what the box binds

In the thesis form (`p.l1.constrain_applied = false`) it binds on μ₁ only. The
thesis says so explicitly, and it has a visible consequence: the realized torque
is `ũ_ff + (L_g̃L_f̃y)⁻¹(μ₁ + μ₂)` and only the μ₁ part was inside the box.
`l1.u_box_excess` reports how far past, per call, rather than leaving it to be
discovered from a plot. With the default the box binds on the total torque (fix
3 above), and the excess is zero by construction.

**The box is sized to the gait, with headroom.** The thesis's 65 Nm is below the
feedforward of every gait in `Results/`, and a μ₁ box that cannot deliver the
feedforward starves the inner QP; the adaptation then chases the tracking
failure instead of the model error and runs away. So `ch4_load_gait` raises
`p.l1.u_max` to **1.25×** the gait's own peak torque (244 Nm on `posture_195`).
The 25% is measured (at ε = 0.35, three steps): at mass scale 1.5 a box exactly at the 195 Nm peak ran
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
| 2 | `ch4_compare_controllers(..., 'robust')` | §4.1.4, Cases I–III |
| 2b | `ch4_compare_controllers(..., 'case4')`, with its own bounds | §4.1.4, Case IV (Fig 4.5) |
| 3 | `ch4_compare_controllers(..., 'l1')` | §4.2.4 |
| 3b | `ch4_load_study` — an unknown torso load, random and fixed | §4.2.4 (Fig 4.11) |
| 4 | `ch4_plot_uncertainty`, `ch4_plot_load` | Figs 4.2–4.11 |
| 5 | `ch4_animate` | — |

Stage 0 stops the run if the gait's own collocation residuals, re-evaluated on
the `M/V/G` currently on the path, exceed `p.verify_tol` (see §6 for why a gait
file cannot be trusted to know which dynamics it was solved on).

Stage 1 is not optional and comes first: a robust controller run outside its own
bound is not a robust controller, it is an aggressive one. Its result is also
what runs. Unless you pass `rclf.delta1_max` / `rclf.delta2_max`, `ch4_main` sets
the bounds to 1.2× the Cases I–III maxima measured along a nominal rollout of
the loaded gait. Case IV (scale 3) is not folded in, since it would size the
robust law for a perturbation Cases I–III never apply. It runs as its own sweep
(stage 2b) with bounds measured the same way at scale 3: Δ₁max = 434.2 and
Δ₂max = 0.800. For a uniform scale ‖Δ₂‖ = |1/s − 1| exactly, so at s = 3 it is
2/3, and the Cases I–III bound of 0.514 would run Case IV outside the law's own
hypothesis. `ch4_compare_controllers` prints a note whenever a robust law is
about to do that.

**The sweeps run 25 steps.** Three steps — the horizon of the thesis's
Figs 4.2–4.4 — separate a controller that falls at once from one that does not,
but not one that converges from one that drifts. At the old ε = 0.35 every
three-step row of the robust law looked converged, and over 25 steps it fell in
step 21 of Case I with a perfect model. 25 is also Fig. 4.6's horizon. The full
study takes about ten minutes at 1 kHz, plus about two for Case IV and four
for the load study, and `n_steps` = 3 gives the thesis-style figures. The GIFs keep their own
`anim_steps` (default 4), because a GIF is a fixed 150 frames.

```matlab
ch4_main('n_steps', 3)                      % thesis-style three-step figures
ch4_main('presets', {'case4'}, 'animate', false)   % Case IV alone
ch4_main('presets', {'load'}, 'animate', false)    % the unknown-load study alone
ch4_main('presets', {'l1'}, 'n_steps', 5)
ch4_main('rclf.delta2_model', 'matrix')     % nested fields take dotted names
ch4_main('rclf.boundary_layer', 0)          % the exact robust law, chatter included
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

The headline claim (Remark 4.6) is not "the robust/adaptive controller is
better on average". It is that the controller's convergence behavior is
**unchanged across perturbations** while the baseline's degrades. So read
**down** each controller across scales, not across controllers within a scale.
`Vend/Vmx` makes that explicit. Remark 4.7 adds a within-scale claim for Case I,
where there is no model error at all. There the robust law should track
*better* than both baselines, because it is defending a worst case that is not
happening. B should do slightly worse than A, because its box sits just below
A's peak.

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
perfect model draws barely more than the feedforward itself (196.4 vs 195.0 Nm
at ε = 0.35, where this was measured over three steps): at 0.8× (157 Nm)
`clfqp_con` fell in step 3 and `rclfqp_con` in step 1, while at the 195 Nm floor
both walk all three. Pass `opts.u_box` to fix the boxes by hand.

**`min Fz` is the column that says whether any of it is walking.** It is the
*true* model's normal force. The stance contact is integrated as a pin, so a
controller can demand that the ground pull the foot down and still complete
steps in simulation; a negative entry means those steps are not realizable.

**With κ = 0, read `Vend/Vmx` loosely for `rclfqp_con`.** The exact law's
worst-case term chatters at the sample rate, so V at the final sample lands
anywhere in an order-of-magnitude band. The same Case I run gave 0.021 at a
195.0 Nm box and 0.19 at 195 Nm plus a floating-point hair. Steps, max‖η‖ and
min Fz survive that. The default boundary layer removes the chatter, and with it
this caveat.

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
7.6e-9 (periodicity) and 1.2e-6 (‖η⁺‖). **25 steps per run**, ε = 0.20, control
at 1 kHz, `Δ₁max = 279.1` and `Δ₂max = 0.514` (measured along a nominal rollout,
×1.2), robust boundary layer κ = 1. L₁ with the three fixes of §4: plant-input
predictor at a = 800, Γ = 1e5, both estimates, rows on the applied torque. Result
set `Results/ch4_{robust,l1}_2026-09-13_20-45-50/`; Case IV, with its own
bounds, `Results/ch4_case4_2026-09-14_13-50-53/`; the load study
`Results/ch4_load_2026-09-14_14-40-48/`.

Each cell reads **steps · max‖η‖ · min Fz** over the whole run. A fall is named
by the step `ch4_simulate` rejected.

### §4.1.4 — the robust CLF-QP holds the orbit where the baselines lose it

| controller | Case I, ×1 (box 195 Nm) | Case II, ×1.5 (box 1098 Nm) | Case III, ×0.7 (box 419 Nm) |
|---|---|---|---|
| A `clfqp` (min-norm) | 25 · 0.25 · 32 N | 25 · 16.6 · **−134 N** | 25 · 14.5 · **−1734 N** |
| B `clfqp_con` (box) | 25 · 0.29 · 50 N | 25 · 16.4 · 73 N | 25 · 15.1 · 1 N |
| **C `rclfqp_con`** | 25 · 1.13 · 50 N | 25 · **4.11** · 90 N | 25 · **4.38** · 23 N |

Nothing falls within 25 steps at ε = 0.20, the baselines included. What changed
is how the baselines fail, not whether they do:

- **Tracking:** under either perturbation A and B reach max‖η‖ 14–17.
- **Timing:** their step times wander — 0.184–0.269 s at ×1.5 and 0.277–0.346 s
  at ×0.7, against a nominal 0.2733.
- **Contact:** A needs the ground to pull its foot down on 12.5% / 3.7% of the
  samples.

The robust law holds max‖η‖ at 4.1 / 4.4 and its step times within
0.264–0.284 s. It keeps the true normal force positive throughout (≥ 90 N and
≥ 23 N), and at ×1.5 it uses at most 511 Nm of a 1098 Nm box. Its per-step CLF
peak starts about two orders of magnitude below both baselines' and is still
4–14× below by step 25 (Figure 1). That is Remark 4.6's contrast. The ×1.5 box is large because the
rule reads A's peak over the run, 1372 Nm over 25 steps.

**Its convergence is not unchanged, though, and Remark 4.7 holds only early.**
The robust law's per-step V peak grows across the sweep's 25 steps:

| case | V peak, step 1 → step 25 |
|---|---|
| ×1.5 | 0.0073 → 0.54 |
| ×0.7 | 0.0031 → 0.44 |
| ×1 (perfect model) | 5e-5 → 0.024 |

In Case I it starts as the tightest tracker. Over steps 1–3 its max‖η‖ is 0.09
against A's 0.25, and its V peaks stay below A's through step 5. By step 10
they are above A's, and the run ends at max‖η‖ 1.13. At the old ε = 0.35 the
same growth dropped the robot in step 21 of Case I. Measured there, the growth
survives the boundary layer, the contact rows, the box and D₂ = 0, and it
vanishes with D₁ = 0, so it comes from the robust term itself. The gait is not the cause:
`posture_195`'s hybrid zero dynamics are stable, δ²_zero = 0.746. At ε = 0.20
the growth is slower, not gone.

**The growth levels off.** Run past the sweep's horizon, with each case's own
bounds and box, the robust law settles within 30–50 steps onto a bounded
plateau. It keeps walking at nominal speed, and no case falls:

| case | steps run | plateau: mean per-step V peak · max‖η‖ | settled by | speed, last 10 steps |
|---|---|---|---|---|
| ×1 (perfect model) | 60 | 0.065 · 1.7 | step ~45 | 1.548 m/s |
| ×0.7 | 60 | 0.63–0.74 · 5.5 | step ~30 | 1.552 m/s |
| ×1.5 | 60 | 0.86 · 5.6 | step ~35 | 1.575 m/s |
| ×3 (Case IV) | 120 | 2.9–3.0 · 10.9 | step ~50 | 1.563 m/s |

So the drift is a transient toward ultimate boundedness, not a divergence. The
25-step maxima of 4.1 and 4.4 above are on their way to 5.6 and 5.5. Two things
follow:

- **Remark 4.6 does not hold in the long run.** The steady error grows with the
  perturbation (1.7 / 5.5 / 5.6 / 10.9), though gracefully. Even the 60-step
  maxima under perturbation stay below what the baselines reach within 25 steps
  (14.5–16.6).
- **The price of robustness is a steady error with a perfect model.** In Case I
  the robust law settles at max‖η‖ 1.7, against A's 0.25 over 25 steps.
  Remark 4.7's advantage is lost by step 10 and does not return within 60
  steps.

**The robust law's three cases walk one orbit.** Figure 4 is the analogue of
Fig. 4.6: 25 steps of each controller in each case. The thesis describes three
different orbits there. Here the robust law's three cases lie on top of each
other, and they have to. A uniform mass scale leaves the hybrid zero dynamics
exactly invariant: `ch4_test_model` asserts both halves, that q̈_drift and the
post-impact velocity are independent of s. So a controller that enforces the
virtual constraints walks the nominal orbit whatever s is.

The baselines' orbits do move, by up to about 6° of torso pitch at ×1.5. That
displacement is their tracking error, not a new orbit of the robot. Only the
torso-load study, which changes M non-uniformly, can move the orbit itself.

### §4.1.4, Case IV — the robust law survives ×3 where both baselines fall, and does not stop

Case IV runs as its own sweep with bounds measured at scale 3, ×1.2:
`Δ₁max = 434.2` and `Δ₂max = 0.800`, against an exact ‖Δ₂‖ of 2/3. The box
follows the same rule: 0.8 × A's 2243 Nm peak = 1794 Nm. The thesis used 300 Nm
on its 32 kg robot. Figures 2 and 3 of the result set are Fig. 4.5.

| controller | Case IV, ×3 (box 1794 Nm) |
|---|---|
| A `clfqp` | **falls in step 4** · 19.6 · 5 N |
| B `clfqp_con` | **falls in step 4** · 20.2 · 142 N |
| **C `rclfqp_con`** | 25 · **1.86** · 96 N |

**A and B fail as the thesis says.** Both take three steps whose durations
collapse (0.26 s, then 0.14 and 0.12–0.13 s) while V climbs to 28–30. Their
fourth step is shorter than the 0.15 m floor.

**C regulates the outputs with a slight degradation, as the thesis says — over
the sweep's 25 steps.**

- **Tracking:** max‖η‖ 1.86, against 1.13 in Case I.
- **Contact and torque:** the true normal force stays ≥ 96 N, the torque peaks
  at 1243 Nm so the box never binds, and no QP fails.
- **Residual error:** output y₂ holds a steady offset of 1.5–2°, and the torso
  cycle settles about 1° lower in pitch than the nominal one (Figure 4, and
  0.46–3.59° over the last five steps against 1.69–4.55°). That is tracking
  error, since a uniform scale cannot move the orbit.

**It does not slow to a stop.** The thesis reports walking that "slows down after
several steps to a complete stop", which it attributes to an unstable orbit. C
averages 1.576 m/s over the first three steps and 1.541 m/s over the last five,
against a nominal 1.563, and a 120-step run is still at 1.563 m/s. The orbit
explanation cannot apply to this model. A uniform mass scale leaves the hybrid
zero dynamics exactly invariant, and `posture_195`'s are stable.

**Past step 25 the degradation stops being slight.** Between steps 25 and 50,
C's max‖η‖ grows to about 10.9, and it holds there to step 120 (the plateau
table above). So at ×3 the robust law walks indefinitely, with a steady error
about twice the ×0.7 and ×1.5 plateaus.

### §4.2.4 — with the three fixes, L₁ holds its convergence where the baseline does not

| controller | ×1 (box 244 Nm) | ×0.7 (box 244 Nm) | ×1.5 (box 366 Nm) |
|---|---|---|---|
| A `clfqp` | 25 · 0.25 · 32 N | 25 · 14.5 · **−1734 N** | 25 · 16.6 · **−134 N** |
| B `l1` | 25 · **0.14** · 50 N | 25 · **4.81** · **−362 N** | 25 · **2.17** · **−210 N** |
| C `l1_con` | 25 · **0.15** · 50 N | 25 · **4.58** · 23 N | 25 · **4.51** · 90 N |
| *C, thesis form* | *25 · 0.38 · 49 N* | *25 · 14.5 · −1403 N* | ***falls in step 14*** · *18.7 · −978 N* |

The last row is §4.2 as written: `predictor = 'thesis'`, Γ = 1e4, the box on μ₁
alone. It comes from `Results/ch4_l1_2026-09-13_18-26-49/`.

**With a perfect model L₁ is at least the CLF-QP.** Its max‖η‖ is 0.14 against
0.25, and its CLF between impacts is about four times lower. The estimate is not
quite zero (‖θ̂‖ ≤ 9), because the 1 kHz hold makes the output acceleration
drift within each period. L₁ cancels part of that intersample error.

**Under either perturbation it holds:**

- **Steps:** both laws walk all 25 steps.
- **Tracking:** max‖η‖ 2.2–4.8, against the baseline's 14.5–16.6.
- **Timing:** step times stay within 0.256–0.291 s (the baseline wanders over
  0.184–0.340 s).
- **Convergence:** the per-step CLF peak does not grow. For `l1_con` at ×1.5 it
  moves 0.078 → 0.17 over the run, where the robust law's grows 0.0073 → 0.54
  (§4.1.4).

In Figure 1 both L₁ curves sit one to two orders of magnitude below the baseline
in both perturbed cases, stationary from the first steps. That is §4.2.4's
claim: convergence unaffected by the perturbation, and much better than the
standard CLF-QP.

**The estimate now measures the uncertainty.** The median ‖θ̂‖/‖θ‖ is 0.95–0.99
in both perturbed cases. In the thesis form θ̂ overshot 3–5×.

**`l1_con` is realizable; `l1` is not quite.** With the rows on the applied
torque, `l1_con` keeps the true normal force at ≥ 23 N (×0.7) and ≥ 90 N
(×1.5). `l1` has no rows, and still asks the ground to pull on 1.6% / 1.3% of
the samples, against 11–14% in the thesis form.

**×1.5 sits near an edge, so the sweep behind the default is reported, not just
its winner.** Cells read steps · max‖η‖ over 25 steps, at ×1.5, with the box
scaled and default projection balls:

| predictor rate a | 500 | 632 (critical) | 700 | **800** | 900 |
|---|---|---|---|---|---|
| `l1` | falls in step 16 | falls in step 16 | 25 · 2.2 | **25 · 2.2** | 25 · 2.3 |
| `l1_con` | 25 · 6.3 | 25 · 12.6 | falls in step 12 | **25 · 4.3** | 25 · 2.8 |

`l1` needs a ≥ 700, and there it walks every time. `l1_con` walks at four of the
five rates. Doubling the projection balls changes which runs fall, and at the
default balls `l1_con`'s α estimate presses its ball (‖α̂‖ 208.5 of 210). So at
×1.5 `l1_con`'s result belongs to this tuning rather than to a margin. Its
max‖η‖ here is measured on every tenth solver point, hence 4.3 against the
table's 4.51.

**Past 25 steps `l1_con` holds on a mass scale, and plain `l1` does not.** The
same runs with the same parameters and boxes, extended to 60 steps (120 for
`l1` at ×1.5). Ranges are over 10-step blocks:

| case | controller | steps | max‖η‖ per block | mean per-step V peak | true contact |
|---|---|---|---|---|---|
| ×1 | `l1` / `l1_con` | 60 / 60 | 0.09–0.14 / 0.09–0.15 | ≈ 3e-4, falling slightly | Fz ≥ 50 N |
| ×0.7 | `l1` | 60 | 3.8–5.5 | 0.15–0.35 | **Fz < 0 on 1.3–2.7% in every block** |
| ×0.7 | `l1_con` | 60 | 4.2–5.4 | 0.25–0.51 | Fz ≥ 23 N |
| ×1.5 | `l1` | **falls in step 76** | 2.0–2.2 to step 40, then 3.3, 4.2, 6.5, 9.4 | 0.15 → 2.2 | Fz < 0 on 0–1.4% |
| ×1.5 | `l1_con` | 60 | 4.5 in steps 1–10, 2.0–2.3 from step 21 | 0.44 → 0.16 | Fz ≥ 90 N |

- **`l1_con` holds its convergence.** It is flat with a perfect model, bounded
  without a trend at ×0.7, and at ×1.5 its 25-step max of 4.5 turns out to be
  the early transient: it settles at max‖η‖ ≈ 2.1 from step 21. The true normal
  force stays positive in every block.
- **Plain `l1` at ×1.5 fails slowly.** It runs as well as `l1_con` for 40
  steps. Then its α̂ reaches the projection bound (206 of 210), and the error
  grows block by block until it falls in step 76. At ×0.7 its demand for a
  negative normal force is persistent, not a transient.
- **Even with a perfect model α̂ creeps**, from about 1.6 to 5–6 over 60 steps,
  absorbing the sample-and-hold error. That is far from its bound, but it does
  not settle.

### §4.2.4, Fig 4.11 — L₁ carries 94% of body weight for 25 steps where a box-matched CLF-QP cannot, but not for 60

The robot carries a mass at the hip that no controller is told about. It is
redrawn every step from 0–70 kg (Fig 4.11a), or fixed at 23 / 35 / 46 kg
(Fig 4.11b). Those are the thesis's 0–30 kg and 10 / 15 / 20 kg scaled by this
robot's mass (×74/32), so each load is the same share of body weight. The load
is a point mass at the torso base, with no rotational inertia of its own.

Two design choices, both measured:

- **A second baseline.** Besides §4.2.4's three controllers, the study runs
  `clfqp_con` with the same box and contact rows as `l1_con`. Under load the
  unconstrained CLF-QP's peak torque is 2–3.5 times `l1_con`'s, so a
  comparison against it alone cannot say whether the adaptation or the torque
  is doing the work.
- **The box.** It is 1.25 × the loaded robot's own feedforward peak along the
  orbit. For a uniform scale this rule is exactly the L₁ preset's. Sized by
  mass ratio instead, `l1_con` fell within three steps at 46 and 70 kg.

**A load is not a small mass scale.** Along the orbit:

| load | 23 kg | 35 kg | 46 kg | 0–70 kg (at 70) |
|---|---|---|---|---|
| ‖Δ₂‖ · its isotropic part | 1.99 · 0.19 | 2.30 · 0.23 | 2.47 · 0.26 | 2.72 · 0.30 |
| min eig(I + Δ₂) | 0.49 | 0.41 | 0.37 | 0.30 |
| loaded feedforward peak → box | 288 → 360 Nm | 330 → 412 Nm | 366 → 458 Nm | 445 → 556 Nm |

A mass scale has Δ₂ = (1/s − 1)I, isotropic and below 1. A hip load distorts
the input gain far more in some output directions than in others. Its
eigenvalues stay real and positive, so no direction reverses, but ‖Δ₂‖ > 1 at
every load, which is outside anything the robust CLF-QP's bound can cover. The
study therefore leaves the robust law out, as the thesis does.

Each cell reads steps · max‖η‖ · min Fz over 25 steps, with the true normal
force computed under the load each step carried:

| controller | 0–70 kg random (box 556 Nm) | 23 kg (box 360 Nm) | 35 kg (box 412 Nm) | 46 kg (box 458 Nm) |
|---|---|---|---|---|
| A `clfqp` (no box) | 25 · 8.95 · 38 N | 25 · 5.79 · 39 N | 25 · 5.61 · 18 N | 25 · 5.33 · **−29 N** |
| A′ `clfqp_con` (box, rows) | **falls in step 10** · 8.45 · 30 N | **falls in step 7** · 20.3 · **−13 N** | **falls in step 3** · 6.92 · 48 N | **falls in step 4** · 23.9 · **−131 N** |
| B `l1` | **falls in step 15** · 8.98 · **−306 N** | 25 · 2.63 · 75 N | 25 · 5.06 · **−164 N** | 25 · 9.23 · **−274 N** |
| **C `l1_con`** | 25 · 9.06 · 6 N | 25 · **2.41** · 102 N | 25 · 4.86 · **−23 N** | 25 · 5.22 · **−51 N** |

**Over 25 steps, at the same torque budget, the adaptation is what walks.**
`l1_con` completes all 25 steps of every case, the random 0–70 kg draw
included, which is the thesis's claim. `clfqp_con`, with the same box and rows,
falls in every case within 3–10 steps. The unconstrained `clfqp` also walks
everything, by drawing 2–3.5× the peak torque.

**Over 60 steps the random load brings `l1_con` down.**

| random 0–70 kg, 60 steps | seed 11 (the draws above) | seed 1 | seed 2 | seed 3 |
|---|---|---|---|---|
| `l1_con` (box 556 Nm) | **falls in step 43** | **falls in step 33** | **falls in step 55** | **falls in step 57** |
| `clfqp` (no box) | 60 | 60 | 60 | 60 |

- **`l1_con` falls under every draw sequence,** and its α̂ sits at the
  projection bound (≥ 207 of 210) in nearly every 10-step block. The falls
  follow no single pattern in the load: the steps before them include both
  drops and rises.
- **The unconstrained CLF-QP walks all 60 steps of all four sequences,** at
  max‖η‖ 3.7–12.2 per block. Over the first 25 steps it drew 2–3 times
  `l1_con`'s torque.
- **Under the fixed 46 kg load `l1_con` holds for 60 steps** (max‖η‖ 4.2–5.2
  per block, true Fz negative on 0.07% of two blocks' samples), but also with
  α̂ at its bound.

So the thesis's load claim survives its own horizon here, not a longer one. An
estimate pinned at its projection bound is also what preceded plain `l1`'s fall
at ×1.5. Larger balls are not the fix: at ten times the size, the random case
fell by step 17. The other sequences are `p.load_seed` 1, 2 and 3; the default,
11, is the sequence the load study draws.

### Why L₁ falls in the long runs, and what does not fix it

**The drift is real but it is a symptom.** θ̂ = α̂‖η‖ + β̂ is redundant: while
‖η‖ stays near a value, every split of θ between α̂ and β̂ predicts equally well,
so nothing holds the split. Over long runs the two estimates drift into
near-opposite directions (median cosine −0.91 at ×0.7 and −0.95 at ×1.5), and
α̂ presses its projection bound.

**The falls are the estimator loop outrunning the sample rate.** Three of the
five long-run L₁ falls were re-simulated sample by sample (plain `l1` at ×1.5 in
step 76, `l1_con` under two random load sequences). Each starts within a few
milliseconds of a footstrike that leaves ‖η‖ at 12–14. The prediction-error loop
runs at about √(Γ + Γ_α‖η‖²) rad/s, which at ‖η‖ = 13 is about 4100 rad/s:
4.1 rad per 1 ms sample, past the RK4 advance's stability limit of about 2.8.
Within 1–3 samples θ̂ reaches about 2700–2900 against a true θ of 30–300, the
adaptive torque swamps the QP, and the step collapses. `ch4_test_l1` check 10
isolates the loop at ‖η‖ = 13.

What was tried, over the nine long runs above: 60 steps of `l1_con` at ×1, ×0.7,
×1.5 and 46 kg, 120 steps of `l1` at ×1.5, and the four random load sequences.

| change | runs that fall, of 9 | what it costs |
|---|---|---|
| none | 5 | — |
| leakage pulling α̂ to zero, 2 / 10 / 50 /s | 6 / 3 / 3 | `l1_con` at ×1.5 falls in steps 14–19 at every rate: at ×1.5 α̂ does real work within a step |
| θ̂ kept continuous across footstrikes / θ̂ folded into β̂ there | 6 / 4 | continuity drops ×0.7 in step 6 |
| cap on α's regressor, 1.5 / 1 / 0.5 rad per sample | 3 / 2 / 1 | at 0.5, `l1_con` at ×1.5 tracks at max‖η‖ 6.0–7.4 per block instead of 2.0–4.5, and the random-load runs that survive pass through excursions to 18–48 |
| 0.5 ms control period, no cap (6 of the 9 runs) | 1 of 6 | plain `l1` at ×1.5 walks 120 flat steps; `l1_con` at ×1.5 takes a transient to 19; twice the runtime |
| normalized adaptation, textbook form (loop held at √Γ, 0.32 rad per sample) | 7 | falls where the law as written walks: ×1.5 in step 37, 46 kg in step 4 |
| normalized adaptation, ceiling 0.5 / 0.75 / 1 / 1.5 / 2 rad per sample | 2 / 0 / 0 / 2 / 2 | at 0.75 and 1, `l1_con` at ×1.5 has a typical step (median per-step max‖η‖) of 2.15 / 2.39 against 1.85 as written (4.44 with the cap at 0.5); at 1 it takes a transient to 16.5 |

**Normalized adaptation** (`p.l1.normalized_rate`, `ch4_l1_deriv`) divides both
adaptation laws by m² = max(1, (Γ + Γ_α‖η‖²)/(κ/Δt)²). That holds the
estimator loop at no more than κ rad per sample, while θ̂ keeps α̂‖η‖ in full
and the law is unchanged wherever the loop was already slower. `ch4_test_l1`
check 11 holds check 10's error and shows it settling at κ = 1 and in the
textbook form. On the robot it works as a window. The textbook form, which holds the loop at √Γ at every error, makes
the loop overdamped under the predictor rate a = 800 and slows the estimate
in ordinary walking: in the runs with uncertainty the median per-step peak of
m² runs 8–53. At 1.5 and 2
normalization rarely engages, and the falls return in runs that fell as
written. The limit that matters in closed loop, about 1 rad per sample, is
well under the RK4 advance's own 2.8; this analysis does not derive it.

**Out of sample.** κ was picked on the nine runs above, so the candidates were
rerun on six further 0–70 kg sequences (`p.load_seed` 4–9, `l1_con`, 60 steps):

| six new random-load sequences | falls | typical step, median per-step max‖η‖ | largest max‖η‖, runs that walk |
|---|---|---|---|
| law as written | 4 (steps 21, 32, 47, 49) | 2.92–3.76 | 7.5–7.8 |
| cap on α's regressor, 0.5 | 3 (steps 24, 30, 44) | 3.35–7.96 | 12.8–39.9 |
| normalized, 0.75 | 1 (step 52) | 3.34–5.53 | 7.9–49.8 |
| normalized, 1 | 2 (steps 11, 44) | 3.28–3.92 | 7.7–26.7 |
| 0.5 ms control period, neither | 0 | 3.11–3.37 | 7.1–16.0 |

- **Across all fifteen long runs, normalization at 0.75 falls once**, against 9
  as written and 4 with the cap. The one fall is seed 6, which the law as
  written walks, and which brings down all three 1 kHz limits tried on it.
  Both kinds of limit slow the estimator most when the error is largest,
  which is also when a changed load has to be learned. That is the likely
  reason; it has not been checked sample by sample.
- **What survives is not always clean.** Two of the new sequences pass
  through excursions to max‖η‖ 20–50 at 0.75 and still walk.
- **Only the faster loop is both safe and flat.** At 0.5 ms the law as written
  fell once in twelve runs and never left max‖η‖ 16 on the new sequences.

The regressor cap and normalization stay in the code as
`p.l1.alpha_regressor_rate` and `p.l1.normalized_rate`, both off by default.
The leakage and footstrike variants were removed after these runs. Every
other result in this section runs the law as written at 1 kHz.

**Beyond the lightest load, L₁ does not track better than the unconstrained
baseline:**

| `clfqp` / `l1_con` | 0–70 kg | 23 kg | 35 kg | 46 kg |
|---|---|---|---|---|
| typical step: median per-step max‖η‖ | 3.10 / 3.00 | 2.38 / 1.97 | 2.29 / 3.96 | 2.47 / 4.07 |
| peak torque [Nm] | 1797 / 556 | 729 / 360 | 1146 / 412 | 1598 / 458 |
| per-step peak ‖u‖, mean [Nm] (Fig 4.11b) | 914 / 458 | 643 / 378 | 1031 / 467 | 1451 / 521 |

- **23 kg:** `l1_con` more than halves the worst error and tracks the typical
  step better, on about half the torque.
- **35 and 46 kg:** its worst step matches the baseline's, but its typical step
  is worse.
- **Random load:** the two are level.

As in the thesis's Fig 4.11b, L₁'s torque grows with the load (378 → 467 →
521 Nm per-step peak).

**Contact.** `l1_con` keeps the true normal force positive throughout, except on
3 and 7 samples of about 15 000 at 35 and 46 kg (minimum −23 and −51 N). The
unconstrained `l1` asks the ground to pull on 0.2% of samples at 35 kg and 11.2%
at 46 kg. Under the random load it falls in step 15, the first step after the
load drops from 60 kg to 12 kg. The estimates do follow the uncertainty: median
‖θ̂‖/‖θ‖ is 1.00–1.05 for `l1_con`.

**The thesis form is worse again.** Measured separately at the same boxes, §4.2's
formulation (thesis predictor, Γ = 1e4, box on μ₁ only) reached max‖η‖ of about
27 at 46 kg and about 103 at 70 kg, with θ̂ overshooting θ 2–4×. Its
constrained law fell in step 3 at 70 kg.

**Loads reshape the gait, but these runs cannot say by how much** (Figure 3).
Under the same load the two walking controllers settle on very different torso
cycles:

- the unconstrained baseline's is 2–3° lower in pitch than nominal: −1.2…2.1°
  over the last five steps, against 1.72…4.55°;
- `l1_con`'s grows to 0.7…6.1°, with pitch rates reaching +65°/s against the
  nominal +25.

A load changes M non-uniformly, so unlike a mass scale it can move the orbit
itself. But the gap between two controllers under the same load is tracking
error, and separating it from the load's own effect would take a controller
that is told the true load.

### Known rough edges

- **Remark 4.4, measured.** The robust law's friction rows hold the *nominal*
  prediction at μ ≤ 0.4. Over 25 steps the true demand is p95 **0.22 / 0.12 /
  0.26** in Cases I–III, with peaks of 0.25 and 0.39 in Cases I–II. At ×0.7 it
  briefly reaches 1.49: a momentary slip the controller cannot see, which is the
  gap Remark 4.4 warns about. With the exact law the ×0.7 case needed p95 1.24
  over just three steps. That was the chatter, not the lighter robot: a uniformly
  scaled robot following the same motion needs the same force *ratio*.
- `QPfail` counts are small — `rclfqp_con` 0 / 0 / 17, `l1` 2 / 1 / 1 and
  `l1_con` 0 / 0 / 0 samples of ~16 700. The slack keeps the QP feasible in principle, so these are
  `quadprog` exit-flag failures, and the fallback saturates the feedforward.
- Numbers quoted elsewhere in this document and in code comments as design
  justifications — the κ sweep, the robust box rule, the L₁ box headroom — were
  measured at ε = 0.35, mostly over three steps, and say so where they appear.
- `ch4_forces` re-solves each controller at every ODE output point, while the
  simulation holds each torque for 1 ms. For a chattering law (κ = 0) those
  torques differ pointwise, but not in what the contact columns say: on the
  torque-only runs of §3, recomputing min Fz under the torques actually held
  matched `ch4_forces` to within 53 N on every run.
- Results in `Results/` stamped between 2026-09-02 and 2026-09-12 23:37 predate
  one or more of this section's fixes — a stale gait, copied boxes, falls
  counted as steps, contact rows stripped — and should not be cited. Robust
  results stamped before 2026-09-13 14:48 ran the exact law and chatter,
  everything stamped before 18:26 that day ran ε = 0.35 over three steps, and L₁
  results stamped before 20:45 ran the thesis-form L₁. The 20-27-10 set is the
  fixed L₁ at a = 632, superseded by the sweep above.

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
    ch4_l1_advance.m        RK4 over one control period, fed what the plant got
    ch4_l1_opts.m           the L1 design options, resolved once (old p -> thesis)
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
    ch4_run_entry.m         one run -> one table row; per-step forces under a changing load
    ch4_compare_controllers.m   the §4.1.4 (Cases I-III, Case IV) / §4.2.4 sweeps; measures the robust boxes
    ch4_load_study.m        Fig 4.11: random and fixed unknown torso loads, loaded-feedforward boxes
    ch4_plot_uncertainty.m  Figs 4.2, 4.3, 4.4, 4.5, 4.6, 4.8, 4.9, 4.10
    ch4_plot_load.m         Figs 4.11a/b and the torso phase portrait under load
    ch4_animate.m           controllers racing on one perturbed robot
    ch4_draw_robot.m        the stick figure, shared by ch4_animate and ch4_plot_load
  Test/
    ch4_test_model.m        true-vs-nominal split, the Δ terms, forces under a per-step load
    ch4_test_rclf.m         the robust guarantee, sampled over the ball
    ch4_test_l1.m           projection, error dynamics, filter, behavior
    ch4_test_all.m
```
