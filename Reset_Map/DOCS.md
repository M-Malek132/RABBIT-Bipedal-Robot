# Reset_Map/

Completes the discrete transition **after** the impact map: swaps the roles of
the stance and swing legs and re-plants the new stance foot on the ground, so
the next step starts from a canonical single-support configuration.

| file | signature | role |
|------|-----------|------|
| `relabel_state.m` | `x_new = relabel_state(x)` | swaps stance joints `(q1,q2)`↔swing joints `(q3,q4)` and the corresponding velocities. Pure index permutation. |
| `rabbit_reset_map.m` | `x_reset = rabbit_reset_map(x_post_impact)` | calls `relabel_state`, then shifts the vertical base coordinate so the **new** stance foot sits exactly on the ground (`pz += P_st(q).z`, since foot height = `−pz − C`). The horizontal shift is intentionally disabled (`px` accumulates so the robot walks forward continuously in the world frame). |

> **Sign note:** the vertical re-plant uses `+ P_st(...)` (matching
> `config/make_initial_state.m`). An earlier `−` doubled the foot-height error;
> the effect is small in practice because the impact fires with the swing foot
> already at `z ≈ 0`.

## Pipeline role

```
rabbit_impact_map(x_end) ──▶ relabel_state ──▶ re-plant ──▶ x_start of next step
```

Applied by every multi-step simulator and by the periodicity constraint:
`x_next = rabbit_reset_map(rabbit_impact_map(x_end))`.
