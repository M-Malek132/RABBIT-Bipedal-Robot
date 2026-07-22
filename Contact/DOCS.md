# Contact/

The **discrete** side of the hybrid model: detecting foot strike and applying
the impact (velocity reset). Position `q` does not change at impact — only
velocities.

| file | signature | role |
|------|-----------|------|
| `rabbit_impact_event.m` | `[value, isterminal, direction] = rabbit_impact_event(t, x)` | ODE **event function** (guard). `value = P_sw(q)` height; triggers (terminal, falling `direction = -1`) when the swing foot reaches the ground. Ignores strikes in the first `t < 1e-3` s so the foot can lift off first. |
| `impact_event_wrapper.m` | `[value, isterminal, direction] = impact_event_wrapper(t, xi)` | thin wrapper so the same event works when the ODE state is the **augmented** `xi = [x; torque_cost]` used by the HZD simulators — it slices `xi(1:2*nq)` and calls `rabbit_impact_event`. |
| `rabbit_impact_map.m` | `[x_plus, impulse] = rabbit_impact_map(x_minus)` | **momentum-conserving impact.** Solves the saddle-point system `[M −J_swᵀ; J_sw 0][dq⁺; λ] = [M·dq⁻; 0]` (post-impact swing foot has zero velocity). Returns `x_plus = [q; dq_plus]`; the optional 2nd output `impulse` = λ is the **contact impulse** `[N·s]` (used by `gait_forces` for the Table 3.1 impact-impulse constraint). |

## Pipeline role

```
swing phase ──(rabbit_impact_event fires)──▶ rabbit_impact_map ──▶ Reset_Map/rabbit_reset_map
```

Used by `Simulation/simulate_hzd_gait*`, `Simulation/simulate_n_steps`,
`Trajectory_Optimization/hzd_constraints` (for the periodicity constraint), and
`Visualization/animate_hzd_result`.
