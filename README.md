# The-RABBIT-Robot

A MATLAB framework for modeling, simulation, control, and trajectory optimization of the RABBIT planar five-link biped robot.

---

# Overview

This project implements a modular hybrid dynamical framework for the RABBIT robot, including:

- Robot modeling and kinematics
- Symbolic rigid-body dynamics
- Constrained contact dynamics
- Hybrid impact/reset maps
- Walking simulation
- Controllers
- Hybrid Zero Dynamics (HZD)
- Trajectory optimization
- Visualization and animation

The architecture is designed for:
- research
- education
- nonlinear control
- hybrid locomotion studies
- trajectory optimization experiments

---

# Robot Model

The RABBIT robot is a planar five-link underactuated biped consisting of:

- stance leg
- swing leg
- torso

Generalized coordinates:
```math
q =
\begin{bmatrix}
p_x & p_z & q_1 & q_2 & q_3 & q_4 & q_5
\end{bmatrix}^T
