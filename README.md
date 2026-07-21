# RABBIT Bipedal Robot — MATLAB Simulation & HZD Gait Framework

![MATLAB](https://img.shields.io/badge/MATLAB-Robotics-orange)
![License](https://img.shields.io/badge/License-MIT-blue)
![Status](https://img.shields.io/badge/Project-Research-green)

A modular MATLAB framework for **modeling, simulation, control, and trajectory
optimization** of the **RABBIT planar five-link biped**. The centerpiece is a
working **Hybrid Zero Dynamics (HZD)** gait-optimization pipeline that produces
dynamically feasible, periodic walking gaits and animates them.

![Rabbit Walking](Results/rabbit_animation.gif)

---

## Robot model

RABBIT is a planar, five-link, **underactuated** biped: a torso plus two legs,
each with a hip and knee. Generalized coordinates (7 DOF):

```
q = [ px, pz, qt, q1, q2, q3, q4 ]'
      └─ floating base ─┘  └ stance ┘ └ swing ┘
```

| symbol   | meaning                                   | actuated? |
|----------|-------------------------------------------|-----------|
| `px, pz` | torso base position (world frame)         | no        |
| `qt`     | torso pitch angle                         | no        |
| `q1, q2` | stance-leg hip, knee                      | **yes**   |
| `q3, q4` | swing-leg hip, knee                       | **yes**   |

Full state `x = [q; dq]` (14×1). Four actuators drive `q1..q4`; the floating
base `(px, pz, qt)` is unactuated → **one degree of underactuation** in single
support. Input map `B = [0₃ₓ₄; I₄]` (see `Dynamics/input_matrix.m`).

**Physical parameters** (baked into the symbolic derivation, *not* a separate
config file — see `Dynamics/rabbit_energy_model_generalized_Lagrange.m`,
nested `Mass_Properties`):

- Torso: 10 kg, length 0.75 m
- Each leg link (thigh, shank): 5 kg, length **0.5 m**
- Gravity: 9.8062 m/s²

> ⚠️ **World-frame Z is UP-positive** (ground = 0, above ground = positive),
> but `pz` is stored **down-positive** in the generalized coordinates. The
> world-frame hip height is therefore `-pz`. Mixing these conventions has been
> the source of several sign bugs; kinematic helpers (`P_st`, `P_sw`, `Tt`,
> `get_body_points`) all return **world-frame, up-positive** heights.

### Continuous dynamics

```
M(q)·q̈ + V([q;q̇]) + G(q) = B·u
```

- `M(q)` — mass/inertia matrix (`Dynamics/M.m`)
- `V([q;q̇])` — Coriolis/centrifugal vector (`Dynamics/V.m`)
- `G(q)` — gravity vector (`Dynamics/G.m`)

During single support the stance foot is pinned (holonomic contact
constraint); `Dynamics/rabbit_constrained_dynamics.m` solves the resulting
KKT system using the stance Jacobian `J_st`.

### Hybrid model

Walking is a hybrid system: **continuous** swing-phase dynamics punctuated by
a **discrete** impact when the swing foot strikes the ground.

- Guard / event: swing-foot height = 0 (`Contact/rabbit_impact_event.m`)
- Impact map: momentum-conserving velocity reset (`Contact/rabbit_impact_map.m`)
- Reset map: relabel stance/swing legs + re-plant on ground (`Reset_Map/rabbit_reset_map.m`)

---

## Dependencies

- **MATLAB** (developed against R2021a-era Symbolic Math Toolbox output)
- **Symbolic Math Toolbox** — only to *regenerate* the dynamics (`M`, `G`, `V`,
  `J_st`, …). The generated `.m` files are committed, so day-to-day use does
  **not** require it.
- **Optimization Toolbox** — `fmincon` (SQP) for gait optimization. **Required**
  for the HZD pipeline.

---

## Quick start — the HZD gait pipeline (this is the working entry point)

```matlab
startup                              % add all folders to the path
cd Trajectory_Optimization
rabbit_hzd_trajectory_optimization   % optimize one gait, inspect, animate
```

This will:

1. Load parameters, bounds, and `fmincon` options (`hzd_problem_setup`).
2. Build an initial gait (hand-tuned `x0` + B-spline `make_coeffs`), optionally
   **warm-started** from a saved result in `Results/`.
3. Run `fmincon` to find a periodic gait (`hzd_cost` + `hzd_constraints`).
4. Print a full diagnostic breakdown (`inspect_solution`).
5. Save `Results/hzd_result_<timestamp>.mat` and animate 10 steps
   (`animate_hzd_result`).

To sweep walking speed and build a **gait library**:

```matlab
cd Trajectory_Optimization
build_gait_library                   % speed continuation around the seed gait
```

> A known-good seed gait (periodic, ~0.355 m/s) is saved at
> `Results/hzd_result_2026-07-20_12-06-36.mat` and is used as the default
> warm-start / continuation seed.

### ⚠️ Broken legacy entry points

The following predate the HZD pipeline and **do not run as-is** — they call
functions that no longer exist. They are kept for reference only:

| file | problem |
|------|---------|
| `main_demo.m` | calls `make_initial_state` with 3 outputs (it returns 1), `simulate_n_steps` with 4 args (it takes 3), and `animate_rabbit_stepping_stones` (does not exist) |
| `Controller/rabbit_controller.m` | calls `desired_gait` / `desired_gait_velocity` (do not exist) |
| `Dynamics/rabbit_dynamics.m` | calls `D_matrix` / `C_vector` / `G_vector` (do not exist; the real functions are `M` / `V` / `G`) |

Use `rabbit_hzd_trajectory_optimization.m` as the entry point instead.

---

## Data flow

```mermaid
flowchart TD
    subgraph GEN["Symbolic generation (offline, Symbolic Math Toolbox)"]
        ELM["rabbit_energy_model_generalized_Lagrange.m<br/>(+ Jacobians.m)"]
    end
    ELM -->|"generates .m"| DYN

    subgraph DYN["Dynamics/ (generated + hand-written)"]
        KIN["Kinematics: Tt, T1-T4, P_st, P_sw"]
        MAT["EoM terms: M, V, G, DM"]
        JAC["Jacobians: J_st, J_sw, Jdotdq_*"]
        CD["rabbit_constrained_dynamics<br/>hzd_closed_loop_ode"]
    end

    subgraph PHASE["Controller/ + Utilities/"]
        TH["theta_of_q / dtheta_dq_of<br/>(phase variable θ = qt+q1+½q2)"]
        GBP["get_body_points / check_ground_validity"]
    end

    subgraph OPT["Trajectory_Optimization/ (fmincon SQP)"]
        SETUP["hzd_problem_setup"]
        COST["hzd_cost"]
        CON["hzd_constraints"]
        MAIN["rabbit_hzd_trajectory_optimization<br/>build_gait_library"]
    end

    subgraph SIM["Simulation/"]
        SHG["simulate_hzd_gait<br/>simulate_hzd_gait_full"]
    end

    subgraph HYB["Contact/ + Reset_Map/"]
        EV["rabbit_impact_event"]
        IM["rabbit_impact_map"]
        RM["rabbit_reset_map / relabel_state"]
    end

    subgraph VIZ["Visualization/"]
        AR["animate_rabbit"]
        AHR["animate_hzd_result"]
    end

    DYN --> CD --> SHG
    PHASE --> CD
    SETUP --> MAIN
    MAIN --> COST --> SHG
    MAIN --> CON  --> SHG
    SHG --> EV
    CON --> IM --> RM
    MAIN -->|"z_opt"| AHR --> AR
    MAIN -->|".mat"| RES[("Results/")]
```

---

## Repository layout

Each folder has a `DOCS.md` with a file-by-file breakdown.

| folder | role | docs |
|--------|------|------|
| `Dynamics/` | symbolic-generated kinematics & equations of motion, constrained dynamics, closed-loop ODE | [Dynamics/DOCS.md](Dynamics/DOCS.md) |
| `Contact/` | impact event + impact (velocity reset) map | [Contact/DOCS.md](Contact/DOCS.md) |
| `Reset_Map/` | post-impact leg relabeling + re-plant | [Reset_Map/DOCS.md](Reset_Map/DOCS.md) |
| `Controller/` | phase variable + (legacy) controllers | [Controller/DOCS.md](Controller/DOCS.md) |
| `Trajectory_Optimization/` | **HZD gait optimization (main pipeline)** | [Trajectory_Optimization/DOCS.md](Trajectory_Optimization/DOCS.md) |
| `Simulation/` | single-step / gait integrators | [Simulation/DOCS.md](Simulation/DOCS.md) |
| `Utilities/` | body points + ground-validity checks | [Utilities/DOCS.md](Utilities/DOCS.md) |
| `Visualization/` | stick-figure animation + GIF export | [Visualization/DOCS.md](Visualization/DOCS.md) |
| `config/` | stepping-stone terrain + B-spline params (legacy) | [config/DOCS.md](config/DOCS.md) |
| `Test/` | smoke test | [Test/DOCS.md](Test/DOCS.md) |
| `Results/` | saved gaits (`.mat`), plots, animation GIF | — |

---

## Hybrid Zero Dynamics — implementation notes

The gait design follows Westervelt et al., *Feedback Control of Dynamic Bipedal
Robot Locomotion* (Ch. 6), and the RES-CLF thesis framing (see `DOCS` per
folder). Key implementation choices, and *why* they matter, are documented at
the point of use — a summary:

- **Phase variable is LINEAR in q.** `theta_of_q` returns
  `θ = qt + q1 + 0.5·q2 = c·q` with `c = [0 0 1 1 0.5 0 0]`
  (`dtheta_dq_of` returns the constant `c`). This is the absolute stance-leg
  angle (Westervelt HH6). An earlier `atan2`-based phase variable **wrapped at
  ±π** and derailed the optimizer; the linear form is provably identical over
  all physical poses but cannot wrap.
- **Virtual constraints.** The four actuated joints `q(4:7)` track a clamped
  cubic **B-spline** `yd(s)` (`bspline_eval`, using `BSpline` / `BSpline_derivative`),
  with phase `s = (θ − θ⁻) / (θ⁺ − θ⁻)` clamped to `[0,1]`. A PD law
  (`Kp=400, Kd=40`) zeroes the output error inside `hzd_closed_loop_ode`.
- **Decision vector** (`unpack_z`): `z = [ coeffs(:) ; x_start ]`, size
  `nu·n_coeffs + 14 = 4·6 + 14 = 38`.
- **Cost** (`hzd_cost`): integral of squared torque **per unit step length**
  (Westervelt eq. 6.43). Normalizing by step length prevents a degenerate
  "step in place".
- **Constraints** (`hzd_constraints`):
  - *Equality:* one-step **periodicity** (13), stance-foot on ground (1),
    stance-foot zero velocity (2), `θ_end = θ_plus` (1), and **NEC1** average
    walking rate `L_step / T_step = v_des` (1).
  - *Inequality:* no ground penetration, mid-step **swing-foot clearance**
    (NIC3), and a **step-length floor** (prevents the degenerate scuff).
- **Solver** (`hzd_problem_setup`): `fmincon` SQP, **central** finite
  differences (the event-triggered `ode45` makes forward FD noisy),
  `ScaleProblem` **off** (so reported feasibility is comparable to the
  constraint tolerance).
- **Warm starting & continuation.** Cold starts from `x0` are unreliable;
  `build_gait_library` warm-starts each speed from the previously converged
  gait and marches `v_des` in small steps.

> **Known open item:** the stance foot is held only by an *acceleration-level*
> contact constraint, so `ode45` allows slow numerical **constraint drift**
> during a step. `inspect_solution` reports `stance-foot drift in step` to
> quantify it. This is independent of the optimization and does not invalidate
> a converged gait.

---

## Reproducibility

- **Generated dynamics files.** `Tt, T1–T4, P_st, P_sw, M, DM, V, G` are
  produced by `Dynamics/rabbit_energy_model_generalized_Lagrange.m`;
  `J_st, J_sw, Jdotdq_st, Jdotdq_sw` by `Dynamics/Jacobians.m`. **If you change
  any mass/length/inertia (in the nested `Mass_Properties`) or the kinematics,
  you must re-run these generator scripts** to regenerate the committed `.m`
  files. The generated files carry a "generated by Symbolic Math Toolbox"
  header — do not hand-edit them.
- **Saved results.** `Results/hzd_result_<timestamp>.mat` each contain
  `z_opt, p, fval, exitflag, report`. The `12-06-36` result is the reference
  periodic gait (≈0.355 m/s) used as the warm-start seed. `build_gait_library`
  writes `Results/gait_library_<timestamp>.mat` (`gait_library` struct array +
  `p`).
- **Nondeterminism.** `config/make_random_initial_state.m` uses `randn` with
  **no seed** — random initial states are not reproducible run-to-run. The HZD
  pipeline instead uses a **hard-coded** `x0`, so its results *are*
  reproducible. Add `rng(<seed>)` before random sampling if you need
  determinism there.
- **Orphaned cache.** `Trajectory_Optimization/Trajectory.mat` is committed but
  **not loaded by any code path**; treat its contents as unverified
  (`whos('-file', 'Trajectory.mat')` to inspect).
- **Simscape model.** `Dynamics/Test_dynamics_Simscape.slx` is a standalone
  dynamics cross-check, not part of the runtime pipeline.

---

## References

E. R. Westervelt, J. W. Grizzle, C. Chevallereau, J. H. Choi, B. Morris,
*Feedback Control of Dynamic Bipedal Robot Locomotion*, CRC Press, 2007.

---

## Author

**Mohammad Malek** — Robotics & Control Research
· [malekmohammad.com](https://malekmohammad.com)
· [github.com/M-Malek132](https://github.com/M-Malek132)

## License

Released under the **MIT License**.
