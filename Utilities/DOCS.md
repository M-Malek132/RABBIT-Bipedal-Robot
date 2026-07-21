# Utilities/

Kinematic helpers shared across simulation, optimization, and validation.

| file | signature | role |
|------|-----------|------|
| `get_body_points.m` | `pts = get_body_points(q)` | struct of key body-point **world positions** `[x; z]` (Z up-positive): `hip, torso_top, stance_knee, swing_knee, stance_foot, swing_foot`. Uses the same transform chain (`Tt, T2, T4, P_st, P_sw`) as the animator. |
| `check_ground_validity.m` | `report = check_ground_validity(q, dq, tol)` | validates a configuration: flags any body point **below** ground (`z < −tol`), checks the stance foot is planted (`|z| ≤ tol`), and — if `dq` given — that the stance foot isn't slipping (`‖J_st·dq‖`). Returns `report.is_valid`, `report.violations`, `report.max_penetration`. |

> **Convention:** heights are world-frame **up-positive**; penetration means a
> point dips below `z = 0`. `max_penetration` is reported as depth below ground
> (negative ⇒ everything is safely above). This matches the sign handling in
> `Simulation/simulate_hzd_gait.m`.

Used by `config/make_*_initial_state.m`, `hzd_constraints`, and
`rabbit_hzd_trajectory_optimization.m` (initial-state sanity check).
