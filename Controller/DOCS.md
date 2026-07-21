# Controller/

The **phase variable** used by the HZD framework, plus a legacy PD controller.

## Phase variable (used by the HZD pipeline)

| file | signature | role |
|------|-----------|------|
| `theta_of_q.m` | `th = theta_of_q(q)` | the gait "clock": absolute stance-leg angle, **linear** in `q`: `θ = qt + q1 + 0.5·q2`. Monotonic over a step; `s = (θ − θ⁻)/(θ⁺ − θ⁻)` is the normalized phase. |
| `dtheta_dq_of.m` | `g = dtheta_dq_of(q)` | gradient `∂θ/∂q`, the **constant** row `c = [0 0 1 1 0.5 0 0]`. Kept in sync with `theta_of_q`. |

> **Why linear, not `atan2`.** These previously computed `θ` from the geometric
> hip→foot angle via `atan2`, which **wraps at ±π** and repeatedly derailed the
> optimizer (huge, spurious `theta_end` violations). Because both leg links are
> length ½, `atan2(rel_x, rel_z)` is *exactly* `qt + q1 + q2/2` for all physical
> poses (`|q2| < π`) — so the linear form is the same function with the branch
> cut removed, and its gradient is a constant (no division by `|hip−foot|²`).
> This matches Westervelt Hypothesis HH6 (phase = an absolute angle, linear in
> `q`), on which the HZD decoupling-matrix analysis depends.

## Legacy controller (not on the HZD path)

| file | signature | status |
|------|-----------|--------|
| `rabbit_controller.m` | `u = rabbit_controller(x)` | **broken** — computed-torque law that calls `desired_gait` / `desired_gait_velocity`, which do **not** exist in the repo. Referenced by the (also broken) `main_demo.m`. Kept for reference. |

The HZD controller itself is **not** a separate file — it is inlined in
`Dynamics/hzd_closed_loop_ode.m` (B-spline virtual constraint + PD tracking).
