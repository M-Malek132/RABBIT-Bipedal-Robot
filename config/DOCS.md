# config/

Initial-state construction and (legacy) stepping-stone / B-spline configuration.
The HZD pipeline uses only `make_random_initial_state` (indirectly) and a
hard-coded `x0`; the rest supports the older stepping-stone demo path.

| file | signature | role |
|------|-----------|------|
| `make_initial_state.m` | `x0 = make_initial_state()` | deterministic starting state: a fixed pose with the stance foot projected onto the ground (`pz += P_st(q).z`) and velocities projected so the stance foot has zero velocity (`(I − pinv(J)·J)·dq`). Returns **one** output (14×1). |
| `make_random_initial_state.m` | `x0 = make_random_initial_state(p)` | rejection-samples a random ground-valid pose far enough from `theta_plus`; projects foot onto ground and stance-foot velocity to zero. **Unseeded `randn` ⇒ nondeterministic.** |
| `config.m` | `varargout = config(action, ...)` | persistent singleton holding stepping-stone terrain and B-spline params. `config('init')` generates random stones + calls `init_bspline_params`; `config('get')` returns the stored struct. Read by `animate_rabbit` for terrain. |
| `init_bspline_params.m` | `ctrl = init_bspline_params()` | legacy B-spline controller config (7 segments, cubic, `qa_start/qa_end`, PD gains). **Distinct** from the HZD pipeline's spline params in `hzd_problem_setup`. |

> **Not the robot model.** There is no `Model/` folder; physical parameters
> (masses, lengths, inertias) live in `Dynamics/rabbit_energy_model_generalized_Lagrange.m`
> (nested `Mass_Properties`), baked into the generated dynamics.

> **Coordinate-shift caveat:** `make_initial_state` / `make_random_initial_state`
> both project foot height with `pz += P_st(q).z`, matching `rabbit_reset_map`.
