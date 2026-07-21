# Simulation/

Numerical integration of the continuous swing phase. Two families: the **HZD**
simulators (used by the optimizer and its animation) and the **generic**
step/multi-step simulators (used by the legacy demo path).

## HZD simulators (closed-loop, event-terminated)

Both integrate `Dynamics/hzd_closed_loop_ode` with `ode45`, terminating on
`Contact/impact_event_wrapper` (swing-foot strike). State is augmented with a
running torque²-cost integrand: `xi = [x; cost]`.

| file | signature | returns |
|------|-----------|---------|
| `simulate_hzd_gait.m` | `[x_end, total_torque_sq, max_penetration, status, swing_clearance, T_step] = simulate_hzd_gait(coeffs, x_start, p)` | one step, **summary only**. `status > 0` iff the impact actually fired and no NaNs. Dense-samples the step to compute worst ground penetration and mid-step swing-foot clearance. Called by `hzd_cost`, `hzd_constraints`, `inspect_solution`. |
| `simulate_hzd_gait_full.m` | `[t_out, x_out] = simulate_hzd_gait_full(coeffs, x_start, p)` | one step, **full trajectory** (`x_out` is N×14). Used for animation and the stance-foot-drift diagnostic. |

## Generic simulators (arbitrary controller — legacy path)

| file | signature | role |
|------|-----------|------|
| `simulate_one_step.m` | `[t_out, x_out] = simulate_one_step(x0, controller)` | integrate `Dynamics/rabbit_ode` for one step, terminating on `rabbit_impact_event`. `controller` is a handle `u = controller(x)` (empty ⇒ zero torque). Warns on unrealistic angles/velocities. |
| `simulate_n_steps.m` | `[t_all, x_all] = simulate_n_steps(x0, nSteps, controller)` | chain `simulate_one_step` across `nSteps`, applying `rabbit_impact_map` + `rabbit_reset_map` between steps. |

> Note the two families use **different** step logic. The HZD path is the one
> exercised by the optimizer; the generic path backs `Test/` and the (broken)
> `main_demo.m`.
