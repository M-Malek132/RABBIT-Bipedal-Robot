# Visualization/

Stick-figure animation of a trajectory and GIF export.

| file | signature | role |
|------|-----------|------|
| `animate_rabbit.m` | `animate_rabbit(x_traj)` | renders a walking trajectory frame-by-frame (torso + both legs, feet markers), tracks the camera, optionally draws stepping-stone terrain from `config('get')`, and writes `Results/rabbit_animation.gif`. Accepts `x_traj` as 14×N or N×14. Aborts on empty/NaN input. Local helper `rabbit_points(q)` computes joint world positions via `P_st, P_sw, Tt, T2, T4`. |
| `animate_hzd_result.m` | `animate_hzd_result(z_opt, p, nSteps)` | stitches `nSteps` HZD steps into one trajectory — for each step calls `simulate_hzd_gait_full`, then `rabbit_impact_map` + `rabbit_reset_map` to seed the next — and hands the result to `animate_rabbit`. This is what `rabbit_hzd_trajectory_optimization.m` calls at the end. |
| `make_doc_figures.m` | `make_doc_figures()` | regenerates the four result figures embedded in `docs/RABBIT_documentation.tex` into `docs/figures/`: gait stick-figure, collocation-vs-closed-loop overlay, phase portrait, and node torques against the Table 3.1 box. Reads the newest `Results/col_result_*.mat`; re-run after the reference gait changes, then recompile the document. |

## Pipeline role

```
z_opt ─▶ animate_hzd_result ─(per step: simulate_hzd_gait_full + impact + reset)─▶ animate_rabbit ─▶ Results/rabbit_animation.gif
```

> `animate_rabbit` reads optional stepping-stone terrain from the `config`
> singleton; if `config('init')` was never called it falls back to flat ground.
