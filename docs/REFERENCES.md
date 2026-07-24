# References — mapped to the codebase

Curated learning resources for the methods this repo implements, each tied to the
files it explains. Ordered low-effort → high-effort within each section. Pair this
with [GUIDE.md](GUIDE.md) (the conceptual tour) and the per-folder `DOCS.md` files
(the file-by-file reference).

---

## 1. Start here — closest match to this repo

### Matthew Kelly, *An Introduction to Trajectory Optimization: How to Do Your Own Direct Collocation* (SIAM Review, 59(4), 2017)

The single best companion to `Trajectory_Optimization/Collocation/`. Builds from a
toy problem up to **Hermite–Simpson collocation for a bipedal walking robot** — the
exact transcription used here. Ships documented MATLAB code.

- Free PDF: https://www2.compute.dtu.dk/courses/02465/_assets/kelly2017.pdf
- Official (DOI 10.1137/16M1062569): https://epubs.siam.org/doi/10.1137/16M1062569

**Read alongside:**
- [`Collocation/col_constraints.m`](../Trajectory_Optimization/Collocation/col_constraints.m) — the Hermite–Simpson defects (Kelly §4) are lines 52–59.
- [`Collocation/col_dynamics.m`](../Trajectory_Optimization/Collocation/col_dynamics.m), [`col_cost.m`](../Trajectory_Optimization/Collocation/col_cost.m), [`col_pack.m`](../Trajectory_Optimization/Collocation/col_pack.m) — the decision-vector transcription.

### Russ Tedrake, MIT 6.832 *Underactuated Robotics*

The definitive course for the concepts underlying the whole repo: underactuation,
limit cycles + **the Poincaré map and its spectral radius** (`rho`), walking, and
trajectory optimization. Free text, interactive notes, and full lecture videos.

- Interactive text: https://underactuated.mit.edu/
- Code/repo: https://github.com/RussTedrake/underactuated
- YouTube channel: https://www.youtube.com/channel/UChfUOAhz7ynELF-s_1LPpWg
- Full Spring 2009 playlist: https://www.youtube.com/playlist?list=PL58F1D0056F04CF8C

**Chapters that map here:**
- *Limit Cycles* / *Poincaré analysis* → [`poincare_stability.m`](../Trajectory_Optimization/poincare_stability.m), [`Collocation/col_poincare.m`](../Trajectory_Optimization/Collocation/col_poincare.m). Best explanation of what `rho` means and why `rho < 1` ⇒ walkable.
- *Walking* → the hybrid model in [`Contact/`](../Contact) + [`Reset_Map/`](../Reset_Map).
- *Trajectory Optimization* → [`rabbit_hzd_trajectory_optimization.m`](../Trajectory_Optimization/rabbit_hzd_trajectory_optimization.m) (shooting) vs. `Collocation/`.

---

## 2. The HZD method — the core framework

### E. R. Westervelt, J. W. Grizzle, C. Chevallereau, J. H. Choi, B. Morris, *Feedback Control of Dynamic Bipedal Robot Locomotion* (CRC Press, 2007)

The repo's **primary reference** — cited throughout the code as "Westervelt eq. 6.43",
"HH6", "NEC1/NEC5". Covers virtual constraints, the phase variable, hybrid zero
dynamics, and the orbital-stability conditions. RABBIT is *the* robot of this book.

**Read Ch. 3 & 6 alongside:**
- Phase variable (HH6) → [`Controller/theta_of_q.m`](../Controller/theta_of_q.m), [`dtheta_dq_of.m`](../Controller/dtheta_dq_of.m).
- Virtual constraints + cost (eq. 6.43) → [`hzd_cost.m`](../Trajectory_Optimization/hzd_cost.m).
- NEC1 (walking rate) / NEC5 (orbital stability) → the constraints in [`hzd_constraints.m`](../Trajectory_Optimization/hzd_constraints.m).
- Impact map → [`Contact/rabbit_impact_map.m`](../Contact/rabbit_impact_map.m); reset/relabel → [`Reset_Map/rabbit_reset_map.m`](../Reset_Map/rabbit_reset_map.m).

Jessy Grizzle's group has RABBIT/MABEL talks on YouTube — search
"Grizzle bipedal locomotion hybrid zero dynamics" and "RABBIT walking robot".

---

## 3. The CLF-QP controller

### A. D. Ames, K. Galloway, K. Sreenath, J. W. Grizzle, *Rapidly Exponentially Stabilizing Control Lyapunov Functions and Hybrid Zero Dynamics* (IEEE Trans. Automatic Control, 59(4), 2014)

The RES-CLF that [`Controller/rabbit_clf_qp_controller.m`](../Controller/rabbit_clf_qp_controller.m)
implements: the CARE block solution `P = [√3 I, I; I, √3 I]`, the ε-scaling that sets
the convergence rate, and the CLF-QP with relaxation `δ`. Maps directly to the code's
"thesis Eqs. 3.12–3.18" comments.

Aaron Ames also has lecture series on YouTube — search
"Aaron Ames control Lyapunov functions bipedal" (Caltech / AMBER Lab).

**Read alongside:** lines 122–193 of `rabbit_clf_qp_controller.m` (RES-CLF assembly +
the QP), and the dispatcher [`Dynamics/hzd_ode_rhs.m`](../Dynamics/hzd_ode_rhs.m) that
selects `'pd'` vs `'clfqp'`.

---

## Suggested reading path

1. Tedrake *Limit Cycles* chapter — understand `rho` / Poincaré maps (with pictures).
2. Kelly tutorial — the collocation solver you can run here.
3. Westervelt Ch. 6 — the HZD framing behind the constraints.
4. Ames RES-CLF paper — the QP controller.

This ordering goes low-effort → high-effort and follows the repo's two axes
(solver: shooting vs. collocation; controller: PD vs. CLF-QP).
