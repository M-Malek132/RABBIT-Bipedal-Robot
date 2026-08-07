# Chapter 6 — Dynamic Walking on Stepping Stones with Control Barrier Functions

Chapter 3 built a QP that asks for **stability**. Chapter 5 showed how to put a
**safety** row in a QP when the constrained quantity is not directly actuated.
Chapter 6 puts them together on RABBIT, with a constraint that is unusual in
two ways: it is a constraint on a **discrete event** (where the swing foot
lands), and it is enforced by a **continuous** condition that the foot never
violates during the step.

---

## 1. The one idea, and what it really is

Everything Chapter 6 constrains is a **position** — a foot placement, a head
height, a step width. Position constraints `g_b(q) ≥ 0` have **relative degree
two**: `ġ_b` contains no control, so `L_g g ≡ 0` and the Section 5.1 reciprocal
CBF row degenerates to `0·u ≥ …`. `ch6_test_qp` measures this on RABBIT rather
than asserting it (`dg/d(q̇) = 0` exactly, over sampled states).

Section 6.1.1's fix is to build a relative-degree-**one** function out of `g_b`:

```
h_CBF(q,q̇) = γ_b g_b(q) + ġ_b(q,q̇) ≥ 0                              (6.1)
```

`h_CBF ≥ 0` keeps `g_b ≥ 0`: at `g_b = 0` with the constraint about to be
violated we would need `ġ_b < 0`, which makes `h_CBF` negative first. So (6.1)
is strictly stronger than what is wanted, and unlike what is wanted it is
enforceable.

### It is the r_b = 2 Exponential CBF, written out by hand

Ask for `h_CBF` to decay no faster than exponentially:

```
γ_b ġ + g̈ ≥ −γ (γ_b g + ġ)
⟺  g̈ ≥ −(γ_b γ) g − (γ_b + γ) ġ  =  −K_b η_b ,   η_b = [g; ġ]
```

which is Definition 5.1 with poles `(γ_b, γ)`. `ch6_cbf_row` therefore builds
`K_b` by calling `ch5_pole_gain([γ_b γ], 2)` rather than typing the products
out, and `ch6_test_barrier` asserts the row equals `g̈ + K_b η_b ≥ 0` to 4e-15.

So Sections 6.1 and 6.2 are **not two constructions**. The only substantive
difference is the class-K function on the right:

| `p.cbf.form` | condition | source |
|---|---|---|
| `reciprocal` | `ḣ_CBF ≥ −γ h_CBF³` | Section 6.1, from `B = 1/h_CBF` |
| `exponential` | `ḣ_CBF ≥ −γ h_CBF` | Sections 6.2, 6.4 — eq. (6.12) |

They **agree exactly at `h_CBF = 0`** (both reduce to `ḣ_CBF ≥ 0`) and differ by
`γ(h − h³)` elsewhere — checked exactly in `ch6_test_barrier`. Which means the
choice of form does not change what happens at the boundary, where a barrier is
doing all of its work. What it changes is behaviour far from the boundary
(the cubic is far more permissive) and at `h < 0` (the cubic demands a faster
recovery). The default is `exponential`.

---

## 2. The geometry, and what each circle actually proves

Every constraint is a function of one 2-D point measured **from the stance
foot**. `ch6_kin` returns it with both terms of `r̈ = L_f²r + L_gL_f r·u`.

### 2.1 Footstep placement (6.9)

| | circle | sense | at `h_f = 0` |
|---|---|---|---|
| `g_ST1` | `O1` at `(−R1, 0)`, radius `R1 + l_max` | inside | `l_s ≤ l_max` |
| `g_ST2` | `O2` at `(l_min/2, −R2)`, radius `√(R2² + (l_min/2)²)` | outside | `l_s ≤ 0` or `l_s ≥ l_min` |

The second is **conditional** — the barrier does not forbid stepping backwards
past the stance foot; the conclusion `l_s ≥ l_min` needs `l_s > 0` as well.
`ch6_test_barrier` states it that way and checks both implications by sampling
the set (401 points, 0 mismatches, 2 boundary points skipped as
floating-point ties).

Note what the construction buys. **Neither row is a constraint at the impact
instant.** Both hold over the whole step, and the impact-time bound is a
corollary of a set the foot never leaves. An equality constraint at touchdown
would be a boundary-value problem needing a re-plan every step; this needs
neither.

**`R2` is the scuffing knob (Remark 6.2).** `O2` passes through `(0,0)` and
`(l_min,0)` and arches above the ground between them, so "outside `O2`" is
simultaneously the lower step-length bound and a foot-clearance requirement.
The arch height is

```
c = √(R2² + a²) − R2 ,  a = l_min/2      ⟹   R2 = (a² − c²)/(2c)
```

so `R2` is not a free constant — pick the clearance and it follows
(`ch6_R2_from_clearance`, verified against the closed form to 1e-17). A short
stone buys less clearance (`c ≤ a`, since the semicircle is the highest arc
through both ground points): a 10 cm step cannot be given a 6 cm arch by this
construction, whatever `R2` is.

**Where the swing foot starts.** At the start of step *k* the swing foot is the
foot the robot just stepped past:

```
l_f(0) = −L_(k−1)  (behind the stance foot)     l_f(T) = +L_k = l_s   (6.13)
```

Measured on the reference gait: `l_f(0) = −0.353 m`, `g_ST1 = 0.953`,
`g_ST2 = 0.315`. This is why "outside `O2`" is a real constraint: the circle's
left edge is a few centimetres behind the stance foot, so the foot must climb
**over** the arch to reach the branch where it may land. It cannot slide along
the ground from one branch to the other.

### 2.2 Overhead obstacles (6.2), (6.5)

`ceiling` is `g_C = h_r − h_H`, affine. `circle` keeps the head outside a disc of
radius `R1o` resting on the obstacle — the same requirement localised, so the
robot walks tall everywhere except under the obstacle. Implemented as separate
primitives rather than approximating the plane by a large circle: a circle of
radius `R` has curvature `1/R` entering `g̈` multiplied by `|ṙ|²`, and at the
3.9 m/s the swing foot reaches, a "flat" 100 m circle still contributes
0.15 m/s² of spurious barrier acceleration.

### 2.3 Step width (6.19)–(6.20)

`O3` tangent to `w = w_max` from below, `O4` tangent to `w = w_min` from above,
both through the initial foot position, both **containment** constraints, with

```
R = (l3² + Δw²) / (2 Δw)
```

Verified in `ch6_test_barrier`: tangency to 6e-17, passage through the start
point to 1e-16, and containment ⟹ the bound on a 120×120 grid with zero
counterexamples.

Two findings here that are worth stating plainly.

> **(a) The thesis's `O4F ≥ R4` cannot bound `w_f` from below.** At the tangency
> abscissa the *outside* condition selects `w_f ≤ w_min` — the wrong side — and
> more generally the complement of a finite disc contains points arbitrarily far
> below the line in every direction. Containment does give the bound, and it is
> what "the same principle can be applied" means, since `O3`'s principle *is*
> containment. Implemented as containment; `p.width.o4_sense = 'outside'`
> selects the literal reading, and the test asserts **both** directions — that
> containment implies `w_f ≥ w_min`, and that the literal form admits 7200
> sampled points with `w_f < w_min`.

> **(b) (6.19) is undefined when the previous step width leaves the new
> window.** `R = (l3² + Δw²)/(2Δw)` divides by the gap between `w0` and the
> boundary, so it needs `w_min < w0 < w_max` — there is no circle tangent to
> `w = w_max` from below through a point already above it. Since Section 6.3.3
> draws `w_d` randomly from `[12, 33] cm` with ±2.5 cm windows, **consecutive
> draws routinely violate this**, and the containment geometry as written cannot
> cover its own numerical study.
>
> What Chapter 6 does instead: build the circle anyway (floored at a *microns*
> division guard, not a design constant), so `g` starts **negative** and the
> barrier is asked to *recover* rather than to *maintain*. That is legitimate for
> an ECBF and not for a reciprocal CBF — Chapter 5 makes exactly this point:
> a polynomial condition on the derivatives of `h` stays well defined at `h < 0`
> while `B = 1/h` does not. `geom.recovering` reports it on every call, and
> `ch6_test_foot3d` exercises it: `g` from −0.042 to +0.005, `w_s = 0.280` for a
> window of `[0.275, 0.325]` the nominal profile (0.233) misses entirely.
>
> A corollary worth noticing: with a ±2.5 cm window, a containment-legal width
> command is always within 2.5 cm of the previous width, so the window **always
> contains** `w0` and a nominal profile that holds the width constant is already
> inside it. Containment-mode Case 2 is therefore a *maintenance* problem, never
> a *reaching* one. The reaching version is the recovery mode, and it is a
> strictly weaker claim.

---

## 3. What is validated on what

| Section | thesis | here |
|---|---|---|
| 6.1 | RABBIT | **RABBIT**, the Chapter-3 reference gait |
| 6.2 | RABBIT, moving stones | **RABBIT** |
| 6.3 | DURUS, 23-DoF humanoid | **not modelled** — split, see below |
| 6.4 | MARLO | **RABBIT** with a RABBIT gait library |

There is no DURUS model in this repository, and inventing one would make every
number in Section 6.3.3 unfalsifiable — a 23-DoF model nobody can check is worse
than no model, because it produces plots that look like results. Section 6.3 is
therefore split three ways:

- **Exact, no robot**: the geometry above (`ch6_test_barrier`).
- **On a surrogate** (`ch6_foot3d`): the swing foot as a double integrator in
  `(l_f, w_f, h_f)` — the subsystem (6.17)–(6.20) are *written in* — carrying all
  four barrier rows at once under an input bound. Cases 1, 2 and 3 of Section
  6.3.3 reproduce there.
- **Not claimed**: that DURUS's dynamics can supply the required foot
  acceleration, or the specific ranges of Section 6.3.3. Those are properties of
  that robot.

---

## 4. The controller

```
min_{u,δ}  ‖μ‖² + p₁δ²                                              (6.26)
s.t.  V̇ + λV ≤ δ                        CLF        — the only row that bends
      ġ_i,CBF + κ(h_i,CBF) ≥ 0           ECBF ×n_b  — hard
      F^v_st ≥ δ_N                       normal force
      |F^h_st| ≤ k_f F^v_st              friction cone
      u_min ≤ u ≤ u_max                  input saturation
```

Decision variable is **`u`, not `μ`**. Chapter 5's QPs use `μ`; here the torque
box is a plain bound in `u` (dense in `μ`) and the contact wrench is
`λ = λ_drift + λ_in u` straight out of the same KKT solve that produced `f` and
`g`, so the friction cone and the minimum normal force are *exact* affine rows.
Same feasible set either way; one is shorter and better conditioned.
`ch6_test_qp` checks that with `n_b = 0` it reproduces `ch3_ctrl_clf_qp`
bit-for-bit (rel. diff 0.00e+00).

**Only the CLF row is slacked.** When barrier and CLF conflict, safety wins and
tracking waits. The consequence is that the QP can be **infeasible**, which is a
statement about the robot's actuators rather than a defect — Remark 6.6 is about
exactly it.

### What the cost pulls towards, and why it is not `min ‖μ‖²`

(6.26) is written with `μ_ref = 0` — the minimum-norm control that certifies the
CLF rate. That reading is implemented (`p.cbf.reference = 'ff'`) and it is **not**
the default. The measurement that decides it, over 10 steps from the reference
fixed point, with a foothold window so wide it can never bind and **the barrier
rows removed entirely**:

| controller | steps | step length |
|---|---|---|
| PD | **10/10** | 0.3537 every step, 110 Nm |
| min-norm CLF-QP (`'ff'`, n_b = 0) | 4/10 | 0.358 → 0.371 → 0.418 → falls, 156 Nm |
| Chapter-6 QP, n_b = 0 | 4/10 | identical to the line above |
| Chapter-6 QP **with** barriers | 4/10 | identical again |

The last two rows are the point. The Chapter-6 QP with *zero* barrier rows
destabilises the gait exactly as the min-norm CLF-QP does, and putting the
barriers back changes nothing — so the failure is **inherited from the cost**.
Every symptom of it (drifting step length, saturated torque, infeasible samples)
looks like a barrier-tuning problem, and no amount of pole tuning touches it.
Chapter 3 had already measured the same thing about its stage-7 law:
*minimum-norm is not the same as well-behaved*.

Setting `μ_ref` to the PD control makes the program a **safety filter**: where no
row is binding it returns PD exactly (`ch6_test_qp` asserts this to 2.9e-12),
and where one is binding it makes the smallest departure from PD that satisfies
it. Only the objective changes — the CLF row, the barrier rows, the contact rows
and the box are untouched, so nothing (6.26) guarantees is weakened. It is also
the closer reading of what Section 6.4 says the controller does: *"tracks the
outputs corresponding to this gait ... while maintaining all above
constraints"*.

With `'pd'`, the slack penalty `p₁ ∈ {1, 10², 10⁴}` makes **no difference at
all** — 10/10 steps at 0.3537 and 109.6 Nm in every case — because the QP
already sits where it wants to be.

**Infeasibility minimises the violation rather than dropping every row.**

```
min_{u,s}  s + w‖u − u_ff‖²    s.t.  A_b u ≤ b_b + s,  s ≥ 0,  u ∈ box
```

Holding the saturated feedforward instead — the obvious choice — discards even
the rows that *were* satisfiable, and measured, that turns a marginal step into
a fall: the run then reports "the robot fell" when the truth is "the QP was a
few Nm short for three samples". The tiny quadratic term is a tie-breaker, not
an objective; without it the solver returns an arbitrary point of the optimal
face and consecutive samples chatter. `qp.feasible` is still false and
`qp.viol` reports how far short it fell.

**Everything is sampled at 1 kHz, including the baselines.** The performance
reason is Chapter 3's (a QP law is only piecewise smooth and stalls an adaptive
solver — measured there at 51 908 RHS evaluations for 0.5 % of a step). The
reason it applies to the *smooth* baselines too is fairness: a continuously
evaluated controller enjoys an advantage no digital implementation has, and
Table 6.1 is a comparison. `ch6_test_sim` verifies the hold converges — halving
the period shrinks the mid-swing error by 2.33×.

---

## 5. The pole problem

This is the part of Chapter 6 that the chapter does not discuss and that decides
everything it can achieve.

At touchdown the barrier is nearly tight **by construction** — the foot is
landing on the edge of a small stone — while it is still moving at metres per
second. So

```
h_CBF = γ_b g + ġ ≥ 0     requires    γ_b ≥ |ġ| / g
```

On the reference gait (strike rate 3.88 m/s, ±2.5 cm stones) that ratio reaches
**≈ 170 rad/s** for `g_ST2`. But the row demands `g̈ ≥ −K_b η_b`, whose dominant
term near `h_CBF = 0` is `−(γ_b + γ)ġ ≈ γ_b|ġ|` — at `γ_b = 170` and
`|ġ| ≈ 2 m/s` that is **340 m/s²** of swing-foot acceleration, which no torque
box on this robot can supply.

**The two requirements are incompatible for this gait.** What tuning buys is
that the infeasible samples end up confined to the last few milliseconds of the
step, *after* the placement is determined: the foot arrives on the stone and the
certificate lapses at the terminal instant. `ch6_report` counts those samples;
`ch6_step_range` separates **landed** from **landed with the guarantee intact**
for exactly this reason.

Two levers follow. Measured one step from the reference fixed point (which lands
at 0.3537 m), sweeping the desired foothold over `0.15:0.025:0.55` m, torque box
300 Nm, PD reference; "landed" = the foot came down inside the window:

**Poles × stone size:**

| poles | ±2.5 cm stone | ±5 cm | ±7.5 cm |
|---|---|---|---|
| 20 | 0.150 – 0.200 | 0.150 – 0.250 | 0.150 – 0.250 |
| **40** | 0.150 – 0.325 | 0.150 – 0.350 | **0.150 – 0.375** |
| 80 | ragged islands | ragged islands | ragged islands |

Below 20 the barrier fights the nominal gait everywhere; above 40 it does
nothing until it demands the impossible and the range breaks into islands.
Bigger stones help monotonically, which is the geometry talking: `g` at
touchdown grows with the distance past `l_min`, and the required pole is
`|ġ|/g`.

**The torque box.** At 150 Nm the achievable range collapses to a couple of
centimetres either side of the nominal; at 300 Nm it is the table above. The
reference gait itself needs 110 Nm, so the box is not what makes the gait
possible — it is what limits how far the barrier can bend it. Quoting a step
range without quoting the box would be quoting nothing.

### The range is asymmetric, and that is physical

The barrier **shortens** the step over a wide range (down to 0.15 m, −58 % of
nominal) and **cannot lengthen** it past about +6 %. That asymmetry is not a
tuning artifact:

- Shortening means putting the foot down early — the barrier only has to stop
  the foot going further, and gravity is already helping.
- Lengthening means holding the swing foot airborne while the *underactuated*
  stance leg keeps rotating. At 1.174 m/s with a 0.301 s step there is very
  little time to buy, and the one degree of underactuation is exactly the
  direction the controller cannot push.

> **An earlier version of this measurement appeared to show lengthening to
> 0.475 m.** It was an artifact: under `p.cbf.reference = 'ff'` the min-norm law
> drifts the step length *upward on its own* (0.358 → 0.371 → 0.418, with no
> barrier at all), so a long "achieved" step was the instability rather than the
> barrier. The numbers above are taken with the PD reference, where the
> no-barrier baseline holds 0.3537 m exactly for ten steps, so any change in
> `l_s` is attributable to the barrier.

This asymmetry is the concrete form of the limitation Section 6.4 exists to
remove: the CBF alone works over a limited range around **one** nominal gait, so
pair it with a library that supplies a gait whose nominal step length is already
near the target and let the barrier handle only the transient.

### The geometric lever, which matters more than the gains

`∂g_ST2/∂h_f = R2/ρ2` is what couples the barrier to the vertical strike rate.
Working the touchdown condition through:

```
γ_b ≥ |ġ| / g  ≈  ( (R2/ρ2)|ḣ_f| − (a/ρ2)l̇_f ) / ( (a/ρ2) ε )
```

for a foot landing `ε` past `l_min`. A *smaller* `R2` — i.e. a *larger*
requested clearance — makes the barrier mostly horizontal at touchdown and drops
the required pole sharply: at `l_min = 0.35 m` and `ε = 2.5 cm`, a 5 cm
clearance needs `γ_b ≳ 170` while a 10 cm clearance needs `γ_b ≳ 12`. The
nominal gait's own mid-swing clearance is 0.19 m, so a 10 cm arch costs it
nothing.

---

## 6. Limits of this implementation

- **The reference gait is fast.** 1.174 m/s, 0.301 s per step, swing foot
  striking at 3.88 m/s. Every number in §5 scales with that strike rate, and a
  slower gait would have substantially more margin. RABBIT's own experiments ran
  slower.
- **The step-range measurement is one step from the periodic orbit**, which is
  the best case: nothing but the stone is asking the controller to deviate. A
  multi-step run over random stones starts each step from wherever the last one
  ended, and its range is narrower.
- **Section 6.4's gait library is RABBIT's, not MARLO's**, so its grid, its
  achievable range and its Table 6.1 percentages are not comparable to the
  thesis's line by line. `ch6_lib_build` records the *achieved* grid in the file
  and `ch6_lib_alpha` reads it back, so a short library produces more
  extrapolation (which `ch6_report` counts, per Remark 6.6) rather than a wrong
  answer.
- **Table 6.1 is run on a reduced grid.** The thesis's 100 × 10 × 7 × 3 is
  21 000 event-terminated integrations with a QP at 1 kHz. `ch6_table61` prints
  the counts it used and the binomial standard error of each cell (≤ 0.5/√n —
  11 points at n = 20), so a reduced run is never mistaken for the full one and
  differences smaller than the error bar are visibly not differences.

---

## 7. Entry points

```matlab
startup
ch6_test_all                                    % the suite (~15 s)
out = ch6_main;                                 % the studies -> Results/ch6_<stamp>/
R   = ch6_step_range(x0, alpha, ch6_params);    % the §6.1.4 range, measured
lib = ch6_lib_build(ch6_params);                % the §6.4 library (slow)
T   = ch6_table61(p, alpha, lib);               % Table 6.1, reduced grid
```

See [Chapter6/README.md](../Chapter6/README.md) for the file-by-file layout.
