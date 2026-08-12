# Chapter 3 — the big picture

One offline solve produces 24 numbers. The robot never sees the solve — only
those numbers, every timestep, for as long as it walks. Everything else in
[`CH3_OPTIMIZATION_FLOW.md`](CH3_OPTIMIZATION_FLOW.md) and
[`CH3_FUNCTION_MAP.md`](CH3_FUNCTION_MAP.md) is detail underneath this handoff.

![Chapter 3 sequence: Engineer, Optimizer and Verifier exchange messages once
offline to produce alpha, a verified periodic gait. The Engineer deploys alpha
to the Controller, which exchanges state and torque with the Robot every
timestep in a loop, until a foot impact ends the step and the loop
repeats.](figures/ch3_sequence.png)

The top half runs once, before deployment. The loop at the bottom runs every
timestep, forever.

**Only `alpha` crosses.** Not the solver, not time `T` — the controller reads
the robot's own leg angle, not a clock, so it never needs to know the solve
happened at all.

For the full outer workflow (seeding, limit staging, continuation) and the
`fmincon` inner loop, see [`CH3_OPTIMIZATION_FLOW.md`](CH3_OPTIMIZATION_FLOW.md).
For the runtime call chain (`ch3_simulate` → `ch3_step` → `ode45` →
`ch3_control` → `ch3_io_lin`), see [`CH3_FUNCTION_MAP.md`](CH3_FUNCTION_MAP.md).
