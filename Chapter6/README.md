# Chapter 6 — Dynamic Walking on Stepping Stones with Control Barrier Functions

Chapter 3 built a QP that asks for **stability**. Chapter 5 showed how to put a
**safety** row in a QP when the constraint is not directly actuated. Chapter 6
puts those together on RABBIT: the constraint is *where the swing foot lands*,
and the guarantee is that it lands on the stone.

```matlab
startup                 % from the repo root
ch6_test_all            % verify every stage
out = ch6_main;         % run the chapter's studies into Results/ch6_<stamp>/
```

---

## The one idea, and what it really is

Everything Chapter 6 constrains is a **position** — a foot placement, a head
height, a step width — so `g_b(q) ≥ 0` has **relative degree two** and the
Section 5.1 reciprocal CBF cannot be written for it (`L_g g ≡ 0`; the QP row
would contain no control at all). Section 6.1.1's fix is to build a
relative-degree-one function out of `g_b`,

```
h_CBF(q,q̇) = γ_b g_b(q) + ġ_b(q,q̇) ≥ 0                              (6.1)
```

and put a barrier on **that**.

**It is the r_b = 2 Exponential CBF of Section 5.2, written out by hand.** Ask
for `h_CBF` to decay no faster than exponentially:

```
γ_b ġ + g̈ ≥ −γ (γ_b g + ġ)
⟺  g̈ ≥ −(γ_b γ) g − (γ_b + γ) ġ  =  −K_b η_b ,   η_b = [g; ġ]
```

which is Definition 5.1 with poles `(γ_b, γ)` — exactly what `ch5_ecbf_gain`
builds. `ch6_test_barrier` asserts the two rows are identical and that
`K_b` equals `ch5_pole_gain([γ_b γ], 2)`. So Sections 6.1 and 6.2 are not two
constructions; the only substantive difference is the class-K function:

| `p.cbf.form` | condition | where the thesis uses it |
|---|---|---|
| `reciprocal` | `ḣ_CBF ≥ −γ h_CBF³` | Section 6.1, from `B = 1/h_CBF` |
| `exponential` | `ḣ_CBF ≥ −γ h_CBF` | Sections 6.2, 6.4 — eq. (6.12) |

They agree exactly at `h_CBF = 0` (both reduce to `ḣ_CBF ≥ 0`) and differ by
`γ(h − h³)` elsewhere. `ch6_test_barrier` checks that identity too.

---

## The geometry

Everything is written in one two-dimensional point measured **from the stance
foot** — `ch6_kin` returns it, its velocity, and the two terms of
`r̈ = L_f²r + L_gL_f r · u`.

### Footstep placement (6.9), Fig. 6.3

| | circle | sense | at `h_f = 0` gives |
|---|---|---|---|
| `g_ST1` | `O1` at `(−R1, 0)`, radius `R1 + l_max` | inside | `l_s ≤ l_max` |
| `g_ST2` | `O2` at `(l_min/2, −R2)`, radius `√(R2² + (l_min/2)²)` | outside | `l_s ≥ l_min` (given `l_s > 0`) |

`ch6_test_barrier` checks both implications by **sampling the set**, not by
trusting the algebra — a barrier can be perfectly enforced and still fail to
deliver what the chapter says it does, and no dynamic test would notice.

**`R2` is the scuffing knob (Remark 6.2).** `O2` passes through `(0,0)` and
`(l_min,0)` and arches above the ground between them, so "outside `O2`" is the
lower step-length bound *and* a swing-foot clearance requirement. The arch
height at mid-chord is `c = √(R2² + a²) − R2`, `a = l_min/2`, so

```
R2 = (a² − c²) / (2c)          ch6_R2_from_clearance
```

— pick the clearance and `R2` follows. A short stone buys less clearance
(`c ≤ a`), which is a real geometric limit rather than a numerical one.

**Where the swing foot actually starts.** At the start of step *k* the swing
foot is the foot the robot just stepped *past*, so

```
l_f(0) = −L_(k−1)      (behind the stance foot)
l_f(T) = +L_k          = the step length at impact             (6.13)
```

Measured on the reference gait: `l_f(0) = −0.353 m`, `g_ST1 = 0.953`,
`g_ST2 = 0.315`. This is why "outside `O2`" is a real constraint — the circle's
left edge sits a few centimetres behind the stance foot, so the foot must climb
**over** the arch to reach the branch where it is allowed to land. It cannot
slide along the ground from one branch to the other.

### Overhead obstacles (6.2), (6.5)

`ceiling`: `g_C = h_r − h_H`. `circle`: the head stays outside a disc of radius
`R1o` resting on the obstacle — the same requirement, localised, so the robot
walks tall everywhere except under the obstacle.

### Step width (6.19)–(6.20), Fig. 6.12

`O3` tangent to `w = w_max` from below, `O4` tangent to `w = w_min` from above,
both through the initial foot position, both **containment** constraints:

```
R = (l3² + Δw²) / (2 Δw)
```

> **The thesis writes the lower bound as `O4F ≥ R4` — the foot *outside* the
> disc — and that cannot bound `w_f` from below.** At the tangency abscissa the
> outside condition selects `w_f ≤ w_min`, the wrong side; more generally the
> complement of a finite disc contains points arbitrarily far below the line in
> every direction. Containment does give the bound, and it is what "the same
> principle can be applied" means, since `O3`'s principle is containment. It is
> implemented as containment. `p.width.o4_sense = 'outside'` selects the literal
> reading, and `ch6_test_barrier` asserts **both**: that the inside form implies
> `w_f ≥ w_min` on sampled points, and that the outside form admits points with
> `w_f < w_min`.

---

## What is validated on what

| Section | subject in the thesis | here |
|---|---|---|
| 6.1 | RABBIT | **RABBIT**, the Chapter-3 reference gait |
| 6.2 | RABBIT, moving stones | **RABBIT** |
| 6.3 | DURUS, 23-DoF 3D humanoid | **not modelled** — see below |
| 6.4 | MARLO | **RABBIT** with a RABBIT gait library |

**There is no DURUS model in this repository**, and inventing one would make
every number in Section 6.3.3 unfalsifiable. Section 6.3 is therefore split:

- **Checked exactly, no robot** (`ch6_test_barrier`): the tangency and
  through-the-initial-point properties of (6.19); that containment implies the
  width bounds; the derivative chain against finite differences.
- **Checked on a surrogate** (`ch6_foot3d`, `ch6_test_foot3d`): that the
  ECBF-CLF-QP carrying all four barrier rows at once drives a swing foot into
  both windows under an input bound, from a nominal profile aimed elsewhere.
  The surrogate is the swing foot as a double integrator in `(l_f, w_f, h_f)`
  — the subsystem (6.17)–(6.20) are *written in*.
- **Not checked, and not claimed**: that DURUS's dynamics can supply the
  required foot acceleration, or the specific ranges quoted in Section 6.3.3.

---

## Layout

| file | role |
|---|---|
| `ch6_params.m` | every knob, starting from `ch3_params` |
| `ch6_main.m` | run the studies |
| `Model/ch6_kin.m` | the tracked point to second order in `u` |
| `Model/ch6_foot3d.m` | the Section 6.3 surrogate |
| `Control/ch6_bar_circle.m`, `ch6_bar_affine.m` | geometry primitives with time partials |
| `Control/ch6_bar_lift.m` | geometry × kinematics → `g, ġ, L_f²g, L_gL_f g` |
| `Control/ch6_bar_stones.m`, `ch6_bar_obstacle.m`, `ch6_bar_width.m` | the chapter's constraint sets |
| `Control/ch6_cbf_row.m` | one barrier → one linear row in `u` |
| `Control/ch6_ctrl_cbf_clf_qp.m` | the controller (6.12)/(6.26) |
| `Control/ch6_admissible.m` | Corollary 5.2 at the start of a step |
| `Simulation/ch6_step.m`, `ch6_simulate.m`, `ch6_terrain.m` | the stepping-stone harness |
| `GaitLibrary/ch6_lib_*.m` | (6.22)–(6.24): build, load, interpolate |
| `Analysis/` | reports, the chapter's figures, `ch6_step_range`, `ch6_table61` |
| `Test/ch6_test_all.m` | the suite |

---

## Design decisions worth knowing

**No kinematic Hessian anywhere.** Writing `g = G(r)` and differentiating twice
puts all the curvature in `d²G/dr²`, a 2×2 matrix of the **barrier**, which is
closed form for every constraint in the chapter. The kinematics enter only
through `r̈`, which the generated `Jdotdq_*` files give exactly. So nothing in
Chapter 6 is finite-differenced, and there is no symbolic Hessian of `P_sw`.

**The tracked point is a difference, `P_sw − P_st`.** It is what Fig. 6.3 draws;
it cancels the stance-foot drift that the acceleration-level contact constraint
allows (~9e-07 m per step); and it is invariant to the re-plant inside `Δ`.

**Barrier rows are hard; only the CLF row slacks.** That is what (6.26) says.
The consequence is that the QP can be infeasible, which is a statement about the
robot's actuators rather than a defect — see below.

**The cost pulls towards the PD law, not towards zero.** (6.26) is written as
`min ‖μ‖²`; that reading is available (`p.cbf.reference = 'ff'`) and is not the
default. Measured over 10 steps with the foothold window so wide it can never
bind **and the barrier rows removed entirely**: PD walks 10/10 at 0.3537 m every
step, while the min-norm CLF-QP manages 4/10 and falls — and the Chapter-6 QP
with zero barriers reproduces the failure exactly. So the instability is in the
*cost*, not the barrier, and every symptom of it looks like a barrier-tuning
problem that no pole retuning fixes. With `μ_ref = μ_PD` the program is a safety
filter: identical to PD where no row binds (asserted to 2.9e-12), smallest
departure from PD where one does. Only the objective changes.

**Infeasibility minimises the violation rather than dropping every row.**
Holding the saturated feedforward discards even the rows that *were*
satisfiable, and measured, that turns a marginal step into a fall — the run then
reports "the robot fell" when the truth is "the QP was a few Nm short for three
samples". `ch6_ctrl_cbf_clf_qp` instead solves
`min s  s.t.  A_b u ≤ b_b + s, u ∈ box` and reports `qp.feasible = false` with
`qp.viol = s`. The guarantee still lapsed and the run still says so.

**Everything is sampled at 1 kHz, including the baselines.** A QP law is only
piecewise smooth, and an adaptive solver stalls on the active-set kinks
(Chapter 3 measured 51 908 RHS evaluations for 0.5 % of a step). The reason it
applies to the *smooth* baselines too is fairness: a continuously evaluated
controller enjoys an advantage no digital implementation has, and Table 6.1 is a
comparison.

---

## The pole problem, measured

This is the part of Chapter 6 that the chapter does not discuss and that decides
everything.

At touchdown the barrier is nearly tight by construction — the foot is landing
on the edge of a small stone — while it is still moving at metres per second. So
`g → small` and `|ġ|` stays large, and

```
h_CBF = γ_b g + ġ ≥ 0   requires   γ_b ≥ |ġ| / g
```

On the reference gait (strike rate 3.88 m/s, ±2.5 cm stones) that ratio reaches
**≈ 170 rad/s** for `g_ST2`. But the ECBF row demands
`g̈ ≥ −K_b η_b`, whose dominant term near `h_CBF = 0` is `−(γ_b + γ) ġ ≈ γ_b|ġ|`
— so a pole large enough to keep the condition satisfied is also a pole that
demands foot accelerations of hundreds of m/s², which no torque box can supply.

**These two requirements are incompatible for this gait.** With the poles
tuned as far as they go, the infeasible samples are confined to the last few
milliseconds of the step, *after* the placement is already determined. The
placement is achieved; the certificate lapses at the terminal instant.

Two things follow, and both are visible in `p.stones` and `p.cbf`:

- The geometry parameter matters more than the gains. `∂g_ST2/∂h_f = R2/ρ2`
  is what couples the barrier to the 3.88 m/s vertical strike rate, so a
  *smaller* `R2` (a *larger* requested clearance) makes the barrier mostly
  horizontal at touchdown and drops the required pole sharply.
- Stone size matters monotonically, for the same reason: `g` at touchdown grows
  with the distance past `l_min`, and the required pole is `|ġ|/g`.

**The achievable range is asymmetric.** The barrier shortens the step down to
0.15 m (−58 % of nominal) and cannot lengthen it past about +6 %. Shortening
means putting the foot down early, which gravity already helps with; lengthening
means holding the swing foot airborne while the *underactuated* stance leg keeps
rotating, and at 1.174 m/s with a 0.301 s step there is very little time to buy.
That limitation is the concrete form of the one Section 6.4 exists to remove.

Run `ch6_step_range` to measure the achievable range for any gait and setting.
It separates **landed** (the foot came down inside the window) from
**feasible** (and the QP never fell back) from **clean** (and no enabled limit
was violated), because those fail in that order and mean different things.
