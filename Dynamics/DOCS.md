# Dynamics/

Robot kinematics and equations of motion. Most files here are **auto-generated**
by the Symbolic Math Toolbox and should **not** be hand-edited — regenerate them
from the two generator scripts instead (see *Reproducibility* in the top-level
README).

## Generator scripts (run offline; require Symbolic Math Toolbox)

| file | role | produces |
|------|------|----------|
| `rabbit_energy_model_generalized_Lagrange.m` | Master Euler–Lagrange derivation. Defines the transform chain, link masses/inertias (nested `Mass_Properties`: torso 10 kg/0.75 m, links 5 kg/0.5 m, g = 9.8062). Script, no inputs. | `Tt, T1, T2, T3, T4, P_st, P_sw, M, DM, V, G` |
| `Jacobians.m` | Symbolically differentiates foot positions to build contact Jacobians. Script, no inputs. | `J_st, J_sw, Jdotdq_st, Jdotdq_sw` |

## Generated kinematics (homogeneous transforms & foot positions)

All take the config `q` (7×1) and return **world-frame** quantities
(Z up-positive).

| file | signature | returns |
|------|-----------|---------|
| `Tt.m`  | `T = Tt(q)`  | torso frame (4×4) |
| `T1.m`–`T4.m` | `T = T1(q)` … | link frames (4×4): stance thigh/shank, swing thigh/shank |
| `P_st.m` | `p = P_st(q)` | stance-foot position `[x; z]` |
| `P_sw.m` | `p = P_sw(q)` | swing-foot position `[x; z]` |

Because both leg links have length ½, the stance hip→foot vector satisfies
`hip − P_st = cos(q2/2)·[sin(qt+q1+q2/2); cos(qt+q1+q2/2)]`, which is why the
phase variable `θ = qt + q1 + q2/2` is exact (see `Controller/theta_of_q.m`).

## Generated equations-of-motion terms

| file | signature | returns |
|------|-----------|---------|
| `M.m`  | `M = M(q)`         | 7×7 mass/inertia matrix |
| `V.m`  | `V = V([q;dq])`    | 7×1 Coriolis/centrifugal vector |
| `G.m`  | `G = G(q)`         | 7×1 gravity vector |
| `DM.m` | `DM = DM([q;dq])`  | time-derivative of `M` (used to form `V`) |
| `J_st.m` | `J = J_st(q)`    | 2×7 stance-foot Jacobian |
| `J_sw.m` | `J = J_sw(q)`    | 2×7 swing-foot Jacobian |
| `Jdotdq_st.m` | `a = Jdotdq_st(q,dq)` | 2×1 stance `J̇·dq` |
| `Jdotdq_sw.m` | `a = Jdotdq_sw(q,dq)` | 2×1 swing `J̇·dq` |

Equation of motion: `M(q)·q̈ + V([q;dq]) + G(q) = B·u`.

## Hand-written dynamics

| file | signature | role |
|------|-----------|------|
| `input_matrix.m` | `B = input_matrix()` | constant 7×4 actuation map `[0₃ₓ₄; I₄]` (actuates `q1..q4`) |
| `rabbit_constrained_dynamics.m` | `[ddq, lambda] = rabbit_constrained_dynamics(q, dq, u[, foot_ref, alpha, beta])` | **single-support forward dynamics.** Solves the KKT system with the stance-foot holonomic constraint (`J_st·q̈ + J̇dq = 0`). Returns accelerations and contact force `lambda = [F_x; F_z]`. Optional `foot_ref` enables **Baumgarte** constraint stabilization (see the closed-loop-ODE note). Called every ODE step. |
| `hzd_control_and_dynamics.m` | `[ddq, tau, lambda] = hzd_control_and_dynamics(q, dq, coeffs, theta_minus, theta_plus, p[, foot_ref])` | **the control law + constrained dynamics for one state**, factored out of the ODE so `tau` (joint torque) and `lambda` (contact force) can be **recomputed from a trajectory after a sim** (see `Trajectory_Optimization/gait_forces.m`) without duplicating the controller. |
| `hzd_closed_loop_ode.m` | `dxi = hzd_closed_loop_ode(t, xi, coeffs, theta_minus, theta_plus, p[, foot_ref])` | **closed-loop swing-phase RHS.** Thin wrapper: calls `hzd_control_and_dynamics` and appends the running torque²-cost integrand. State `xi = [x; torque_cost]`. Optional `foot_ref` (world stance-foot pin) enables Baumgarte in the dynamics. |
| `clf_qp_closed_loop_ode.m` | `dxi = clf_qp_closed_loop_ode(t, xi, coeffs, theta_minus, theta_plus, p)` | **CLF-QP analogue of `hzd_closed_loop_ode`.** Torque comes from `Controller/rabbit_clf_qp_controller.m` (input-output linearization + control-Lyapunov-function QP) instead of the fixed PD law; returns the QP's constrained-dynamics accelerations and appends the running torque²-cost integrand. Same `xi = [x; torque_cost]` state, so the existing simulate/plot machinery applies. |
| `hzd_ode_rhs.m` | `dxi = hzd_ode_rhs(t, xi, coeffs, theta_minus, theta_plus, p[, foot_ref])` | **Controller dispatcher** — the single place that selects the swing-phase feedback law. `p.controller = 'pd'` (default) → `hzd_closed_loop_ode`; `'clfqp'` → `clf_qp_closed_loop_ode`. `simulate_hzd_gait`/`_full` call this, so the optimizer runs whichever law `p.controller` names. Unset/unknown ⇒ `'pd'` (backward compatible with saved `p`). |
| `rabbit_ode.m` | `dx = rabbit_ode(t, x, controller)` | generic single-support RHS for an arbitrary `controller(x)` (empty ⇒ zero torque). Used by `Simulation/simulate_one_step.m`. |

## Broken / orphaned

- `rabbit_dynamics.m` — *unconstrained* forward dynamics, but calls
  `D_matrix` / `C_vector` / `G_vector`, which **do not exist** (the real terms
  are `M` / `V` / `G`). Not on any working path; only named in `startup.m`'s
  dependency check.

## Non-code

- `Test_dynamics_Simscape.slx` / `.slxc` — standalone Simscape model for
  cross-checking the derived dynamics. Not part of the runtime pipeline.
