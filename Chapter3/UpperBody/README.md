# Chapter3/UpperBody — upper-body wrench on the RABBIT torso

Adds a wrench `W = [Fx; Fz; M]` acting on the torso at a "shoulder" point
(`l_sh` from the hip along the torso axis), entering the EOM as

    M(q) q̈ + V + G = B u + J_st' λ + Q_ub,     Q_ub = J_sh(q)' [Fx; Fz] + e_qt M

`J_sh` has non-zero columns only for `px, pz, qt`, so the wrench acts directly on
x, z, q_t and reaches the legs only through `M(q)` and the stance constraint.
With no `p.ub` field every Chapter-3 result is unchanged (test T1).

| file | role |
|------|------|
| `ch3_ub_defaults.m` | attach `p.ub` (modes `none` / `profile` / `bezier`, bounds, weight) |
| `ch3_ub_wrench.m`   | `Q_ub`, `W`, shoulder point and Jacobian (complex step through `Tt.m`) |
| `ch3_ub_info.m`     | the one place that reads `p.ub`'s switches |
| `ch3_ub_cost.m`     | effort penalty on an optimized wrench |
| `ch3_ub_test.m`     | 7 checks (identity when off, Jacobian, virtual power, Newton balance, packing) |
| `ch3_ub_demo.m`     | A: sensitivity on the orbit, B: closed loop, C: re-optimization |

Patched existing files (all no-ops without `p.ub`):
`Model/ch3_control_affine.m` (Q_ub in the drift column),
`HZD/ch3_zd_point.m` (Q_ub next to G in the zero dynamics),
`Optimization/ch3_col_pack.m`, `ch3_col_unpack.m`, `ch3_col_effective_params.m`,
`ch3_col_bounds.m` (β appended last in z), `ch3_col_cost.m` (+ wrench effort),
`ch3_col_constraints.m` (optional W(0) = W(1)).

Not modelled: impulsive wrenches at impact (the reset map is unchanged, correct
for any bounded W), and the wrench's own dynamics — an internal upper body
(arms) should be an extra link, not a free W.

**Crutch-axis cone** (`'k_cone'`, off by default): `|Fx| <= k_cone * Fz` on every
Bezier coefficient, which guarantees it for the whole curve. Appended as extra
rows after row 19 of `ch3_col_constraints`, only when set.
