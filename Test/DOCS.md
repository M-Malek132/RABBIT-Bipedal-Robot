# Test/

| file | role |
|------|------|
| `test_simulate_steps.m` | Smoke test for the **generic** simulation path: `startup` → `config('init')` → `make_initial_state` → `simulate_n_steps(x0, 10, [])` → `animate_rabbit`. Script, no inputs/outputs. |

> This exercises the legacy generic simulator (zero-torque, `[]` controller),
> **not** the HZD pipeline. It does not assert numeric results — it runs the
> chain and animates, so "passing" means "runs without error and looks right".
>
> There is currently no automated test for the HZD gait optimization; use
> `Trajectory_Optimization/inspect_solution.m` on a saved `z_opt` for that.
