# Visualization/

Stick-figure animation of a trajectory and GIF export.

| file | signature | role |
|------|-----------|------|
| `animate_rabbit.m` | `animate_rabbit(x_traj, opts)` | renders a walking trajectory frame-by-frame (torso + both legs, feet markers), tracks the camera, optionally draws stepping-stone terrain from `config('get')`, and writes a GIF. Accepts `x_traj` as 14×N or N×14. Aborts on empty/NaN input. `opts` is optional — `.gif_file` (default `Results/rabbit_animation.gif`), `.skip` (5), `.delay` (0.03), `.title`, `.name` — and the defaults reproduce the original one-argument behaviour, so older callers are unaffected. Local helper `rabbit_points(q)` computes joint world positions via `P_st, P_sw, Tt, T2, T4`. |
| `animate_hzd_result.m` | `animate_hzd_result(z_opt, p, nSteps)` | stitches `nSteps` HZD steps into one trajectory — for each step calls `simulate_hzd_gait_full`, then `rabbit_impact_map` + `rabbit_reset_map` to seed the next — and hands the result to `animate_rabbit`. This is what `rabbit_hzd_trajectory_optimization.m` calls at the end. |
| `animate_col_result.m` | `animate_col_result(col_file, nSteps, opts)` | animates a **direct-collocation** gait from `Results/col_result_*.mat` (newest by name if omitted) and writes `Results/collocation_animation.gif`. A collocation solution is N node states plus a step time, not a time series, so this rebuilds the trajectory the transcription implies — the cubic Hermite polynomial through `(x_k, f_k)`, whose midpoint equals the `xmid` that `col_constraints` uses — and tiles the periodic step `nSteps` times, advancing `px` by one step length each time (no re-simulation, no accumulated error). Prints step length/speed, the playback factor, and the periodicity defect; warns if that defect exceeds 1e-3. `opts.n_sub` defaults to real-time playback. |
| `make_doc_figures.m` | `make_doc_figures()` | regenerates the four result figures embedded in `docs/RABBIT_documentation.tex` into `docs/figures/`: gait stick-figure, collocation-vs-closed-loop overlay, phase portrait, and node torques against the Table 3.1 box. Reads the newest `Results/col_result_*.mat`; re-run after the reference gait changes, then recompile the document. |

## Pipeline role

```
z_opt ─▶ animate_hzd_result ─(per step: simulate_hzd_gait_full + impact + reset)─▶ animate_rabbit ─▶ Results/rabbit_animation.gif

col_result_*.mat ─▶ animate_col_result ─(Hermite interpolant + periodic tiling)─▶ animate_rabbit ─▶ Results/collocation_animation.gif
```

The two animators answer different questions. `animate_hzd_result` shows what the
**closed loop** does (it integrates the controller, so it diverges when the orbit
is not attracting — which this seed orbit is not, ρ ≈ 3). `animate_col_result`
shows what the **optimizer** produced, open loop; use `col_crosscheck` to compare
the two.

> `animate_rabbit` reads optional stepping-stone terrain from the `config`
> singleton; if `config('init')` was never called it falls back to flat ground.

> **Between-node drift is real, not a rendering artifact.** `col_constraints`
> pins the stance foot only *at* the nodes (to ~1e-12 there). On the current
> N = 20 mesh, h = 55 ms, the Hermite interpolant between nodes lifts the stance
> foot up to ~6 mm and slides it ~2.4 cm, so the animation shows a foot that
> creeps. That is a property of the mesh, and raising N is the fix — do not
> "smooth" it away in the animator.
