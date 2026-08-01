# Understanding RABBIT — A Guided Tour of the HZD Gait Pipeline

*Onboarding for someone who knows robotics/control but is new to **this** codebase.*

This guide teaches you how the code **thinks**, not just what each file does. The
[README](../README.md) and the per-folder `DOCS.md` files are the reference
manual; this is the tutorial that makes them make sense. Read this once,
top-to-bottom, then keep the `DOCS.md` files open as you work.

> **One-sentence summary.** The repo answers a single question — *"what periodic
> walking motion can this 5-link biped execute, and is it any good?"* — by turning
> gait design into a constrained optimization problem (`fmincon`) over the shape
> of a few splines, and then simulating and animating the result.

---

## 1. The 30-second mental model

You are looking at a **Hybrid Zero Dynamics (HZD)** gait-design toolbox for one
specific robot: the planar RABBIT biped. Everything in the repo exists to support
one loop:

```
        pick spline shapes  ─────────────►  simulate one step  ─────────►  score it
        (decision variables)                (ode45 + impact)              (cost + constraints)
              ▲                                                                  │
              └──────────────  fmincon adjusts the shapes to improve the score ◄─┘
```

When `fmincon` converges, you have a set of spline shapes that produce a
**periodic** step: the robot ends the step in the same configuration it started
(after you swap its legs), so it can repeat forever. That converged gait is saved
to `Results/` and animated.

Three ideas do 90% of the conceptual work. If you understand these, the rest is
plumbing:

1. **The robot is a hybrid system** — smooth swing dynamics + an instantaneous
   impact when the foot hits the ground.
2. **Virtual constraints** — instead of commanding joint *torques* directly, we
   command the four motors to *track spline curves*, and we let a PD (or CLF-QP)
   law figure out the torques. The splines' shapes are what we optimize.
3. **The phase variable** — the splines are indexed not by *time* but by a
   "gait clock" `θ` derived from the robot's own pose. This is what makes the
   whole thing self-clocking and time-invariant.

The next three sections unpack exactly these.

---

## 2. Core concept 1 — The robot as a hybrid system

### The coordinates

RABBIT has **7 generalized coordinates** (planar, 5 links):

```
q = [ px, y, qt,  q1, q2,  q3, q4 ]'
      └ floating base ┘ └stance┘ └swing┘
```

- `px, y` — torso base position, `qt` — torso pitch. **Unactuated** (no motor
  attaches the torso to the world).
- `q1,q2` — stance hip & knee, `q3,q4` — swing hip & knee. **Actuated** (4 motors).

Full state is `x = [q; dq]`, a **14×1** vector. Four actuators, seven DOF →
**one degree of underactuation**. That single unactuated DOF is *the whole point*:
a fully actuated robot is a trivial trajectory-tracking problem; the underactuated
torso is what makes walking hard and interesting, and it is what the "zero
dynamics" in HZD refers to.

> ⚠️ **The sign trap that has bitten everyone.** World-frame Z is **up-positive**
> (ground = 0). But `y` is stored **down-positive**. So world hip height is
> `−y`. Every kinematic helper (`P_st`, `P_sw`, `Tt`, `get_body_points`) already
> returns **world-frame up-positive** heights — trust them, and never hand-roll a
> height from `y`.

### The two-phase cycle

Walking is not one dynamical system; it is two, stitched together:

```
   ┌──────────────────────────────────────────────┐
   │  CONTINUOUS  (single support / swing phase)   │
   │  Stance foot pinned to ground; swing foot     │
   │  arcs forward. Governed by an ODE.            │
   └───────────────┬──────────────────────────────┘
                   │  guard:  swing-foot height = 0
                   ▼
   ┌──────────────────────────────────────────────┐
   │  DISCRETE  (impact + reset)                   │
   │  1. impact map:  velocities jump (q is fixed) │
   │  2. reset map:   relabel legs, re-plant foot  │
   └───────────────┬──────────────────────────────┘
                   │
                   └────────► start of next step
```

**Continuous side** (`Dynamics/`, `Simulation/`). During single support the
stance foot is a pin joint. The equation of motion

```
M(q)·q̈ + V([q;q̇]) + G(q) = B·u        (+ the pin constraint J_st·q̈ + J̇·q̇ = 0)
```

is solved as a **KKT system** in `rabbit_constrained_dynamics.m` — accelerations
and the contact force `λ` fall out of one linear solve `A\b`. `M`, `V`, `G`,
`J_st` are all **symbolically generated** files (more on that in §6).

**Discrete side** (`Contact/`, `Reset_Map/`). Two steps, applied in this exact
order (`rabbit_reset_map(rabbit_impact_map(x_end))`):

1. **Impact map** (`rabbit_impact_map.m`) — the swing foot hits and sticks.
   Configuration `q` is **continuous** (positions don't teleport); only
   **velocities** jump, via another KKT solve that conserves momentum while
   forcing the new contact point's velocity to zero.
2. **Reset map** (`rabbit_reset_map.m`) — the foot that just landed *is* the new
   stance foot, so we **relabel** legs (swap `q1,q2 ↔ q3,q4`) and shift `y` so
   the new stance foot sits exactly on the ground. `px` is deliberately **not**
   reset — it accumulates so the robot walks forward in the world.

**Why order matters:** the impact map uses the *pre-relabel* swing Jacobian
`J_sw`. Relabel first and you'd apply the impact to the wrong foot. The
docstring in `rabbit_impact_map.m` shouts this for a reason.

---

## 3. Core concept 2 — Virtual constraints (the thing being optimized)

Here is the conceptual leap that trips up newcomers. **We do not optimize a
torque trajectory.** We optimize the *shape of curves that the joints are asked to
follow*, and a feedback law manufactures whatever torque is needed to follow them.

Concretely, the four actuated joints `q(4:7)` are commanded to track a **desired
profile** `yd(s)`:

```
output error:   y  = q(4:7) − yd(s)        ← we want this driven to zero
```

Each of the 4 desired profiles `yd_i(s)` is a **clamped cubic B-spline** with
`n_coeffs = 6` control points. So the design knobs are `4 × 6 = 24` spline
control-point values, plus the `14`-element initial state — **38 numbers total**.
That vector is `z`:

```matlab
% unpack_z.m — this is the whole "language" fmincon speaks
z = [ coeffs(:) ; x_start ]        % 4*6 + 14 = 38
[coeffs, x_start] = unpack_z(z, p);
```

**How the torque is produced.** Inside the swing-phase ODE
(`hzd_control_and_dynamics.m`), a fixed-gain PD law zeroes the output error:

```matlab
tau = -p.Kp.*e - p.Kd.*de;         % Kp=400, Kd=40 ; e = y − yd, de = ẏ − ẏd
```

That `tau` goes straight into the constrained dynamics. So the optimizer never
sees torques as variables — torque is a *consequence* of the spline shapes and the
tracking law. (There is a fancier alternative — a CLF-QP that solves a small
quadratic program for `tau` every ODE step and can *guarantee* torque limits — see
§8. The default is PD.)

**Why splines, why clamped, why B-splines specifically?**

- **Few parameters, smooth output.** 6 control points give a smooth curve with
  enough freedom to shape a stride, but few enough that `fmincon` over 38 vars is
  tractable.
- **Clamped** knot vector means the **first control point is exactly the curve's
  value at `s=0`**. `make_coeffs.m` exploits this to guarantee `yd(0)` equals the
  initial joint angles exactly — no discontinuity at step start.
- The Cox–de Boor evaluation lives in `BSpline.m` / `BSpline_derivative.m`;
  `bspline_eval.m` wraps them and also central-differences a second derivative
  when the CLF-QP needs it.

---

## 4. Core concept 3 — The phase variable `θ` (the gait clock)

The splines above are indexed by `s ∈ [0,1]`, the **phase**. Where does `s` come
from? Not from a stopwatch — from the robot's own configuration:

```matlab
% theta_of_q.m
theta = qt + q1 + 0.5*q2         %  = c·q,   c = [0 0 1 1 0.5 0 0]
```

Then `s = (θ − θ⁻) / (θ⁺ − θ⁻)`, clamped to `[0,1]`, where `θ⁻` is `θ` at the
start of the step and `θ⁺ = p.theta_plus = 0.3` is the target end. So the gait is
**self-clocking**: the robot's own leg angle tells it how far through the step it
is. This is what makes the controller **time-invariant** — push the robot and it
re-indexes rather than falling off a pre-planned timeline.

`θ` is the **absolute stance-leg angle**, and it is deliberately **linear in `q`**.
This is the single most instructive "war story" in the codebase:

> **Why linear and not `atan2`?** `θ` was once computed as the geometric hip→foot
> angle via `atan2(rel_x, rel_z)`. Because both leg links have length ½, that
> `atan2` is *exactly* `qt + q1 + q2/2` for every physical pose — but `atan2`
> **wraps at ±π**. When the optimizer pushed a leg past that, `θ` jumped from
> `+2.6` to `−3.04` and `fmincon` chased a phantom discontinuity. The linear form
> is the *same function with the branch cut removed*, and its gradient is the
> **constant** row `c` (see `dtheta_dq_of.m`) instead of a nonlinear expression
> that divides by `|hip−foot|²`. This matches Westervelt's Hypothesis HH6, on
> which the entire HZD decoupling analysis depends.

**Takeaway for a newcomer:** `theta_of_q` and `dtheta_dq_of` must always agree
(`dtheta_dq_of` returns `c`, the gradient of `theta_of_q`). If you ever change one,
change the other.

---

## 5. The optimization problem, assembled

Now the three concepts combine into the `fmincon` problem. Every piece lives in
`Trajectory_Optimization/`.

### The decision vector, cost, and constraints

| piece | file | in one line |
|---|---|---|
| **variables** `z` (38×1) | `unpack_z.m` | spline control points + initial state |
| **cost** `J` | `hzd_cost.m` | ∫‖u‖² **per unit step length** |
| **constraints** `[c, ceq]` | `hzd_constraints.m` | periodicity + physics + realizability |
| **setup** `p, lb, ub, options` | `hzd_problem_setup.m` | single source of truth for everything |

**The cost is normalized by step length**, and that `1/L_step` is not cosmetic:

```matlab
% hzd_cost.m
J = total_torque_sq / L_eff;       % L_eff = max(L_step, min_step_length)
```

Without the division, a short step is cheap, so the optimizer is *rewarded for
shrinking the stride* — the "step in place" degeneracy the comments warn about.
Normalizing makes the objective "energy per distance," which is what you actually
want. (There is also a smooth quadratic penalty below a floor, deliberately used
instead of a hard cliff so a gradient always exists — read the comment block in
`hzd_cost.m`; it is a small masterclass in why cost functions must stay
continuous for gradient-based solvers.)

**The constraints** are where the physics lives. `hzd_constraints.m` simulates one
step and then assembles:

*Equalities (`ceq`, 18)* — these define "a valid periodic step":
- **Periodicity (13):** `x_next(2:14) == x_start(2:14)` after impact+reset. This
  is *the* HZD condition — the step returns to its start (excluding `px`, which
  advances).
- Stance foot on ground (1), stance foot zero velocity (2).
- `θ_end == θ_plus` (1) — the step actually reaches the target phase.
- **NEC1 walking rate (1):** `L_step / T_step == v_des` — pins the speed.

*Inequalities (`c`, 8)* — "physically sane and realizable":
- No ground penetration; mid-step swing-foot **clearance** (so the foot doesn't
  scuff early); a **step-length floor**.
- **Orbital stability** `rho < 1` (gated — see §7).
- Four **physical-realizability** limits (torque / friction / min-GRF / impact
  impulse), each individually toggleable.

> **Design pattern to notice:** *gated constraints held at `−1`.* A disabled
> inequality is set to `−1` (trivially satisfied) rather than removed, so `c`
> always has length 8. This keeps `fmincon`'s problem dimensions constant while
> letting you switch constraints on one at a time. This is the mechanism behind
> the whole "enable them incrementally" workflow.

### Solver choices worth understanding

`hzd_problem_setup.m` sets `fmincon` to **SQP** with **central** finite
differences and `ScaleProblem` **off**. Each choice has a reason:

- **Central** FD because the cost/constraints come from an **event-terminated
  `ode45`** — the impact time is a nonsmooth function of `z`, so forward
  differences are noisy. Central differences average out some of that noise.
- **Scaling off** so the reported `Feasibility` number is directly comparable to
  `ConstraintTolerance` (with scaling on, a failed *first* evaluation can poison
  the internal scaling — hence the "keep the failure sentinel the same order of
  magnitude as real constraints" comment in `hzd_constraints.m`, value `10` not
  `1e3`).

---

## 6. Where the dynamics come from (and why you rarely touch them)

Open `Dynamics/M.m` and you'll find a wall of machine-generated trig. **Do not
hand-edit these.** They are *output*, produced offline by two generator scripts
that require the Symbolic Math Toolbox:

```
rabbit_energy_model_generalized_Lagrange.m   ──►  Tt, T1–T4, P_st, P_sw, M, DM, V, G
Jacobians.m                                  ──►  J_st, J_sw, Jdotdq_st, Jdotdq_sw
```

The **physical parameters** (torso 10 kg / 0.75 m, each link 5 kg / 0.5 m,
g = 9.8062) are **baked into the symbolic derivation** — there is no runtime
config file for the robot model. If you change a mass or length, you must **re-run
the generators** to regenerate the committed `.m` files. This is the single most
important reproducibility rule in the repo.

The generated files are committed precisely so that **day-to-day use does not need
the Symbolic Math Toolbox** — only the Optimization Toolbox (`fmincon`).

The hand-written glue that ties generated terms into the pipeline:

- `rabbit_constrained_dynamics.m` — the single-support KKT solve (§2).
- `hzd_control_and_dynamics.m` — control law + dynamics for one state, factored
  out so torques/forces can be **recomputed after a sim** without duplicating the
  controller (this is how `gait_forces.m` measures Table-3.1 quantities).
- `hzd_closed_loop_ode.m` / `clf_qp_closed_loop_ode.m` — the two swing-phase RHS
  variants (PD and CLF-QP).
- `hzd_ode_rhs.m` — the **dispatcher** that picks between them based on
  `p.controller`. This one file is why flipping `p.controller = 'clfqp'` changes
  the *entire* optimization with no other edits.

---

## 7. Periodic ≠ stable — the most important subtlety in the repo

A newcomer's natural assumption is "if the periodicity constraint is satisfied,
the robot walks." **This is false, and understanding why is the key insight.**

The periodicity equality makes `x_start` a **fixed point** of the step-to-step map
`P: x_start → reset(impact(simulate(x_start)))`. But a fixed point can be
**unstable**. The reference seed gait has `max|periodicity| ≈ 0.007` (a
near-perfect fixed point) yet a Poincaré spectral radius `rho ≈ 3`. Simulate it
open-loop and the velocities grow `6 → 6 → 19 → 891 → …` — **it falls over within
a few steps.**

This is Westervelt's **NEC5**: a periodic orbit must also be *attracting*
(`rho < 1`). It's the gait-scale analogue of the "step in place" degeneracy — a
missing constraint the optimizer will happily exploit.

- `poincare_stability.m` builds the 13×13 Poincaré-map Jacobian by
  central-differencing the step map (~26 step sims) and returns its spectral
  radius `rho`. `rho < 1` ⇒ walkable; `rho ≥ 1` ⇒ periodic-but-falls.
- Setting `p.enforce_stability = true` adds the inequality `rho ≤ rho_max` (< 1).

**Why it's gated off by default:** each evaluation is ~26 simulations, and
`fmincon` finite-differences it → hundreds-to-thousands of sims per iteration.
The intended workflow is: **solve a periodic gait fast (stability off), then a
short warm-started final phase with stability on** to drive `rho` below 1.

This same "periodic but unstable" story recurs in the collocation solver
(§8), where two feasible gaits measured `rho ≈ 9.8` and `19.8`. **Feasibility and
stability are different properties, and only stability makes a gait walk.**

---

## 8. The two solvers, and the CLF-QP option

There are **two independent ways** to solve the same gait problem. A newcomer
should know both exist and which to reach for.

**Single-shooting** (`rabbit_hzd_trajectory_optimization.m`, the default entry
point). Integrates a full step with `ode45` inside every `fmincon` evaluation; the
initial state is a decision variable. Simple, directly simulates what will run —
but each evaluation is an ODE solve, and it inherits the stance-foot **drift**
issue (the pin is only enforced at the acceleration level, so it slowly slides
~0.1–0.3 m per step).

**Direct collocation** (`Trajectory_Optimization/Collocation/`, the thesis
method). No ODE in the loop: the dynamics become **Hermite–Simpson defect
constraints** between nodes, the stance foot is pinned at *every* node (killing the
drift), and the virtual constraints are enforced at every node so the gait lies
*on* the zero-dynamics manifold by construction. Faster per-evaluation, more
robust — but more variables and its own bookkeeping (`col_pack`/`col_unpack`).

**The controller axis is orthogonal to the solver axis.** `p.controller` selects
the swing-phase feedback law used *everywhere* (sim, cost, constraints,
`gait_forces`) because they all route through `hzd_ode_rhs`:

- `'pd'` (default) — fixed-gain virtual-constraint PD. Fast.
- `'clfqp'` — an input-output-linearizing **control-Lyapunov-function QP**
  (`rabbit_clf_qp_controller.m`) that solves a small QP for the torque at every
  ODE step, enforcing a hard torque box (and optionally friction) *inside* the
  feedback. A gait optimized this way is torque-feasible **by construction** — but
  a QP-per-ODE-step is slow, so **warm-start it from a PD gait; never cold-start.**

---

## 9. End-to-end trace — one run of the default pipeline

Let's follow `rabbit_hzd_trajectory_optimization` from start to finish. This is the
"single run through the whole pipeline" that makes the pieces click.

```matlab
startup                               % add paths, verify deps, print "Startup complete"
cd Trajectory_Optimization
rabbit_hzd_trajectory_optimization    % the whole thing
```

**Step 0 — Setup.** `hzd_problem_setup()` returns `p`, bounds `lb/ub`, and
`fmincon` options. The script pins `p.v_des = 0.355` (the seed gait's natural
speed, so it's already near-feasible) and enables all four physical constraints.

**Step 1 — Initial guess.** A hand-tuned `x0` (14×1) defines a natural
post-impact pose (swing foot behind stance, 6 cm clearance, `θ = −0.15` so phase
sweeps `−0.15 → 0.30`). `make_coeffs(x0, p)` builds the initial splines by
interpolating from `x0`'s actuated joints to a **relabeled** touchdown pose — the
relabel is a smart guess at what a *periodic* touchdown should look like. Because
the spline is clamped, `yd(0)` equals `x0`'s joints exactly.

**Step 2 — Warm start (recommended).** The script loads a known-good result
(`Results/hzd_result_2026-07-20_12-06-36.mat`) and uses its `z_opt` as `z0`. Cold
starts from `x0` "wandered for 100+ iterations"; warm starting from a periodic gait
is what makes this reliable. It keeps the **current** `p` (so today's constraint
settings apply) and only borrows the vector.

**Step 3 — Solve.**

```matlab
[z_opt, fval, exitflag] = fmincon(@(z) hzd_cost(z,p), z0, [],[],[],[], lb,ub, ...
                                  @(z) hzd_constraints(z,p), options);
```

Every `fmincon` evaluation does the following (this is the heart of the whole
repo):

```
z ──unpack_z──► (coeffs, x_start)
                     │
                     ▼
        simulate_hzd_gait(coeffs, x_start, p)
                     │   ode45( hzd_ode_rhs )  from θ⁻ to swing-foot strike
                     │      each RHS eval:  s = phase(q) → yd(s) via spline
                     │                      tau = −Kp·e − Kd·de
                     │                      [ddq,λ] = rabbit_constrained_dynamics
                     │      terminate on impact_event_wrapper (swing foot z = 0)
                     ▼
             x_end, T_step, penetration, clearance
                     │
        ┌────────────┴─────────────┐
        ▼                          ▼
   hzd_cost                  hzd_constraints
   ∫‖u‖²/L_step         x_next = reset(impact(x_end))
                        ceq = x_next − x_start  (periodicity) + foot + θ_end + NEC1
                        c   = penetration, clearance, step-floor,
                              [stability], torque/friction/grf/impulse
```

`fmincon` nudges `z`, re-evaluates, and repeats until the constraints are
satisfied (periodic step) and the cost is locally minimal.

**Step 4 — Inspect.** `inspect_solution(z_opt, p)` prints the full report: θ
sweep, step length/duration, **walking speed**, cost, every constraint residual,
worst violation, **stance-foot drift**, and **`rho`** (stability). This is the
first thing you read to judge a result.

**Step 5 — Save & animate.** The result is saved to
`Results/hzd_result_<timestamp>.mat` (`z_opt, p, fval, exitflag, report`), then
`animate_hzd_result(z_opt, p, 10)` stitches 10 steps — each is
`simulate_hzd_gait_full` followed by `impact + reset` to seed the next — and hands
the trajectory to `animate_rabbit`, which writes `Results/rabbit_animation.gif`.

**To sweep speed:** `build_gait_library` marches `p.v_des` outward from the seed,
warm-starting each solve from the previous converged gait (continuation).

---

## 10. Common pitfalls (read before you touch anything)

1. **Running the broken legacy entry points.** `main_demo.m`,
   `Controller/rabbit_controller.m`, and `Dynamics/rabbit_dynamics.m` predate the
   HZD pipeline and **call functions that no longer exist**. They are reference
   corpses. The working entry point is
   `Trajectory_Optimization/rabbit_hzd_trajectory_optimization.m`. (`startup.m`
   deliberately does *not* dependency-check the legacy files, because `which`
   would falsely report them `[OK]`.)

2. **Mixing the two Z conventions.** World Z up-positive vs `y` down-positive
   (§2). Always get heights from the kinematic helpers, never from `y` directly.

3. **Assuming periodic means stable (§7).** The most common conceptual error. A
   converged gait with tiny periodicity residual can still have `rho ≫ 1` and fall
   over. Always check `rho` in `inspect_solution`'s output.

4. **Cold-starting.** Both the CLF-QP path and speed continuation *depend* on
   warm-starting from an existing gait. Cold starts are a gamble; from PD, they're
   often intractable for CLF-QP.

5. **Enabling physical constraints all at once, in the wrong order.** The
   reference gait has vertical GRF going **negative** (the foot would lift off),
   which also makes friction `|Fx/Fz|` blow up at the `Fz=0` crossing. Enable
   **GRF first** (get `Fz > 0`), *then* friction, torque, impulse — warm-starting
   each phase. Enabling friction first makes `fmincon` chase a divide-by-zero and
   diverge. (Documented at length in `hzd_constraints`'s call site.)

6. **Copying thesis numbers.** The thesis limits are for **ATRIAS** (63 kg, 50:1
   gears — its `|u| ≤ 5 Nm` is *motor* torque = 250 Nm at the joint). RABBIT is
   ~30 kg with **no gearing** (`u` *is* joint torque). **Measure first**
   (`inspect_solution` prints all four Table-3.1 quantities), then set limits.

7. **Hand-editing generated dynamics.** `M.m`, `G.m`, etc. are Symbolic-Toolbox
   output. Edit the **generator** and regenerate; never touch the generated files.

8. **The `ddq` clamp is load-bearing.** `hzd_control_and_dynamics.m` clamps
   accelerations to `±1000`. The saved gait library was optimized *against* the
   clamped dynamics — removing the clamp breaks its periodicity (`7.5e-3 → 0.35`).
   Don't remove it unless you're re-optimizing from scratch (`p.clamp_ddq = false`).

9. **Nondeterminism in random initial states.** `make_random_initial_state.m` uses
   unseeded `randn`. The HZD pipeline avoids this by using a hard-coded `x0`, so it
   *is* reproducible; add `rng(seed)` before random sampling if you need
   determinism elsewhere.

10. **Trusting `Trajectory.mat`.** It's committed but loaded by **no code path**.
    Treat it as unverified.

---

## 11. Quick file-finder (where do I look for…?)

| I want to… | look at |
|---|---|
| run the working pipeline | `Trajectory_Optimization/rabbit_hzd_trajectory_optimization.m` |
| change a limit, gain, tolerance, or toggle | `Trajectory_Optimization/hzd_problem_setup.m` |
| understand the objective | `hzd_cost.m` |
| understand the constraints | `hzd_constraints.m` |
| see the control law | `Dynamics/hzd_control_and_dynamics.m` |
| see the single-support dynamics | `Dynamics/rabbit_constrained_dynamics.m` |
| understand the impact/reset | `Contact/rabbit_impact_map.m`, `Reset_Map/rabbit_reset_map.m` |
| understand the phase clock | `Controller/theta_of_q.m` |
| judge a saved result | `Trajectory_Optimization/inspect_solution.m` |
| check stability | `Trajectory_Optimization/poincare_stability.m` |
| use the collocation solver | `Trajectory_Optimization/Collocation/rabbit_hzd_collocation.m` |
| change the robot's mass/length | `Dynamics/rabbit_energy_model_generalized_Lagrange.m`, then regenerate |
| animate a result | `Visualization/animate_hzd_result.m` |

Each folder's `DOCS.md` is the authoritative file-by-file reference.

---

## 12. Glossary

| term | meaning |
|---|---|
| **RABBIT** | The planar 5-link underactuated biped this repo models (torso + 2 legs, each hip+knee). |
| **HZD** | *Hybrid Zero Dynamics.* Gait-design method: enforce virtual constraints so the controlled robot evolves on a low-dimensional "zero dynamics" manifold; design the periodic orbit there. (Westervelt et al., Ch. 6.) |
| **Virtual constraint** | A commanded relationship `q(4:7) = yd(s)` that a feedback law drives to hold. The *outputs* `y = q(4:7) − yd(s)` are zeroed. |
| **Phase variable `θ`** | The "gait clock": `θ = qt + q1 + 0.5·q2`, an absolute stance-leg angle, **linear in `q`**. Normalized to `s = (θ−θ⁻)/(θ⁺−θ⁻) ∈ [0,1]`. |
| **Underactuation** | Fewer actuators (4) than DOF (7). The unactuated torso is the reason walking is nontrivial and the source of the "zero dynamics." |
| **Floating base** | `(px, y, qt)` — the torso's unactuated position/orientation in the world. |
| **Single support / swing phase** | The continuous phase: one foot pinned, the other swinging. |
| **Impact map** | Instantaneous, momentum-conserving velocity reset when the swing foot strikes (`q` continuous, `dq` jumps). |
| **Reset map** | Post-impact relabeling of stance/swing legs + vertical re-plant. |
| **Guard / event** | Swing-foot height = 0; terminates the swing-phase ODE (`rabbit_impact_event.m`). |
| **KKT system** | The saddle-point linear system solved for accelerations + contact force under a holonomic constraint (`[M −Jᵀ; J 0]`). |
| **Poincaré map `P`** | Step-to-step map `x_start → reset(impact(simulate(x_start)))`. Its fixed points are periodic gaits. |
| **`rho` (spectral radius)** | Largest eigenvalue magnitude of `P`'s Jacobian. `rho < 1` ⇒ **stable/walkable**; `rho ≥ 1` ⇒ falls. |
| **Periodic vs. stable** | Periodic = fixed point of `P`; stable = *attracting* fixed point (`rho < 1`). Both are required to walk. |
| **NEC1** | Westervelt "necessary condition": average walking rate `L_step/T_step = v_des`. Pins speed. |
| **NEC5** | Orbital stability requirement (`rho < 1`). |
| **NIC3** | Mid-step swing-foot clearance (swing foot stays above ground through the step's middle). |
| **Table 3.1** | The thesis's four physical-realizability limits: torque, friction cone, min vertical GRF, impact impulse. |
| **GRF** | Ground Reaction Force. `Fz > 0` means the foot is in compression (loaded); `Fz < 0` means it's being pulled down, i.e. lifting off. |
| **CLF-QP** | Control-Lyapunov-Function Quadratic Program: an alternative torque law that guarantees CLF descent + hard torque limits, solved as a QP each ODE step. |
| **Single shooting** | Solver that integrates a full step with `ode45` inside each optimizer evaluation (default). |
| **Direct collocation** | Solver that enforces dynamics as Hermite–Simpson defect constraints between nodes; no ODE in the loop. |
| **Baumgarte stabilization** | Adds `2α·ċ + β²·c` to a constraint so violations decay instead of drifting. Present but **disabled** (`α=β=0`) in this pipeline. |
| **Warm start / continuation** | Seeding a solve from a previously converged gait; marching a parameter (e.g. `v_des`) in small steps, warm-starting each. |
| **Decision vector `z`** | `[coeffs(:); x_start]`, 38×1 — what `fmincon` optimizes. |
| **Gait library** | A struct array of converged gaits across a range of speeds (`build_gait_library.m`). |

---

*Primary reference:* E. R. Westervelt, J. W. Grizzle, C. Chevallereau, J. H. Choi,
B. Morris, *Feedback Control of Dynamic Bipedal Robot Locomotion*, CRC Press,
2007 (esp. Ch. 3 & 6).
