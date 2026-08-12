# Chapter 3 — the big picture

One offline solve produces 24 numbers. The robot never sees the solve — only
those numbers, every timestep, for as long as it walks. Everything else in
[`CH3_OPTIMIZATION_FLOW.md`](CH3_OPTIMIZATION_FLOW.md) and
[`CH3_FUNCTION_MAP.md`](CH3_FUNCTION_MAP.md) is detail underneath this handoff.

```mermaid
sequenceDiagram
    participant Eng as Engineer
    participant Opt as Optimizer (offline)
    participant Ver as Verifier
    participant Ctrl as Controller (online)
    participant Robot as Robot dynamics

    Eng->>Opt: hand-tuned seed pose
    Opt->>Opt: search for alpha (Bezier coeffs)
    Opt->>Ver: does this trajectory hold up?
    Ver-->>Opt: no - refine mesh, try again
    Ver-->>Eng: yes - periodic, stable gait

    Eng->>Ctrl: deploy alpha
    loop every timestep
        Robot->>Ctrl: current state x
        Ctrl->>Robot: torque u(x, alpha)
    end
    Robot->>Robot: step ends - foot impact, repeat
```

The top half runs once, before deployment. The loop at the bottom runs every
timestep, forever.

**Only `alpha` crosses.** Not the solver, not time `T` — the controller reads
the robot's own leg angle, not a clock, so it never needs to know the solve
happened at all.

For the full outer workflow (seeding, limit staging, continuation) and the
`fmincon` inner loop, see [`CH3_OPTIMIZATION_FLOW.md`](CH3_OPTIMIZATION_FLOW.md).
For the runtime call chain (`ch3_simulate` → `ch3_step` → `ode45` →
`ch3_control` → `ch3_io_lin`), see [`CH3_FUNCTION_MAP.md`](CH3_FUNCTION_MAP.md).
