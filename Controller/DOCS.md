# Controller/

The **phase variable** used by the HZD framework, the **CLF-QP controller**,
plus a legacy PD controller.

## CLF-QP controller (thesis Chapter 3)

| file | signature | role |
|------|-----------|------|
| `rabbit_clf_qp_controller.m` | `[u, info] = rabbit_clf_qp_controller(q, dq, coeffs, theta_minus, theta_plus, p)` | Input-output-linearizing **control-Lyapunov-function quadratic-program** controller. Drop-in alternative to the fixed-gain PD virtual-constraint law: same HZD outputs `y = q(4:7) − yd(s)`, same acceleration-level constrained dynamics, but the torque solves a QP that enforces the RES-CLF descent `V̇ ≤ −(γ/ε)V` while respecting a hard torque box (and an optional friction cone). Returns the torque `u` and an `info` struct (`ddq`, `lambda`, `V`, `Vdot`, `delta`, decoupling matrix `Adec`, `Lf2y`, `exitflag`). |

**How it works.** Because `h0(q)=q(4:7)` is linear and the phase `θ` is linear in
`q`, the outputs have relative degree 2 and
`ÿ = L_f²y + (L_gL_fy)·u`, with the decoupling matrix `L_gL_fy = Jy·(∂q̈/∂u)`
and `L_f²y = Jy·q̈₀ − κ²·(d²yd/ds²)·(c·dq)²` (`Jy = H − κ·(dyd/ds)·c`,
`κ = 1/(θ⁺−θ⁻)`, `c = dtheta_dq_of`). The accelerations come from the same KKT
system as `rabbit_constrained_dynamics` (Baumgarte off), written as the affine
map `q̈ = q̈₀ + Mᵤ·u`. In transverse coordinates `η=[y; ẏ]` the CARE
`F'P+PF−PGG'P+I=0` has the closed form `P=[√3·I, I; I, √3·I]`; the ε-scaled
`P_ε` gives the rapidly-exponentially-stabilizing CLF `V=η'P_εη`, and the QP
minimizes `½‖u‖² + p_relax·δ²` subject to the (softened) CLF inequality, the
torque box, and optional friction. Tunables live on `p`
(`clf_eps`, `clf_R`, `clf_relax`, `clf_enforce_torque`, `clf_enforce_friction`);
sensible defaults are baked in.

The swing-phase closed loop under this law is
`Dynamics/clf_qp_closed_loop_ode.m` (the CLF-QP analogue of
`hzd_closed_loop_ode`). Validate with `Test/test_clf_qp.m`, which finite-
difference-checks the analytic Lie derivatives and the CLF descent against a
saved gait.

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
