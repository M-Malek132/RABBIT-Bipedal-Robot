RABBIT Bipedal Robot – MATLAB Simulation Framework








A modular MATLAB framework for modeling, simulation, control, and trajectory optimization of the RABBIT planar five‑link biped robot.

This project implements hybrid dynamical walking and modern nonlinear control tools commonly used in locomotion research.

Demo
Example simulation of the RABBIT robot performing a walking motion.



Overview
This repository provides a structured implementation of the RABBIT robot including:

Robot modeling and kinematics
Symbolic rigid‑body dynamics
Constrained contact dynamics
Hybrid impact/reset maps
Multi‑step walking simulation
Controllers (PD, Feedback Linearization, Hybrid Zero Dynamics)
Hybrid Zero Dynamics (HZD) framework
Trajectory optimization using direct collocation
Visualization and animation tools
The architecture is designed for:

robotics research
nonlinear control studies
hybrid dynamical systems
trajectory optimization experiments
robotics education
Robot Model
The RABBIT robot is a planar five‑link underactuated biped robot consisting of:

torso
stance leg (hip and knee)
swing leg (hip and knee)
The generalized coordinates are

q = [px, pz, qt, q1, q2, q3, q4]ᵀ

where

px, pz – base position
qt – torso angle
q1, q2 – stance leg joints
q3, q4 – swing leg joints
The robot dynamics follow a hybrid model that includes continuous dynamics during swing phase and discrete reset maps at foot impact.

Repository Structure
text
RABBIT-Bipedal-Robot
│
├── main_demo.m
├── startup.m
├── README.md
├── rabbit_animation.gif
│
├── Model/
│   ├── parameters.m
│   ├── packParameters.m
│   └── rabbit_kinematics.m
│
├── Dynamics/
│   ├── D_matrix.m
│   ├── C_vector.m
│   ├── G_vector.m
│   ├── rabbit_dynamics.m
│   └── rabbit_constrained_dynamics.m
│
├── Contact/
│   ├── foot_positions.m
│   ├── rabbit_impact_event.m
│   └── rabbit_impact_map.m
│
├── Reset_Map/
│   └── rabbit_reset_map.m
│
├── Controller/
│   ├── rabbit_controller.m
│   ├── rabbit_pd_controller.m
│   ├── rabbit_fl_controller.m
│   ├── rabbit_hzd_controller.m
│   └── rabbit_virtual_constraints.m
│
├── Simulation/
│   ├── hybrid_simulation.m
│   ├── simulate_one_step.m
│   └── simulate_n_steps.m
│
├── Trajectory_Optimization/
│   ├── direct_collocation.m
│   ├── optimize_gait.m
│   ├── gait_constraints.m
│   └── rabbit_hzd_trajectory_optimization.m
│
├── Visualization/
│   └── animate_rabbit.m
│
├── Utilities/
│   ├── finite_difference.m
│   ├── saturation.m
│   └── save_animation_gif.m
│
├── Test/
│   ├── test_kinematics.m
│   ├── test_jacobians.m
│   └── test_simulate_one_step.m
│
└── Results/
Quick Start
Clone the repository:

text
git clone https://github.com/M-Malek132/RABBIT-Bipedal-Robot.git
cd RABBIT-Bipedal-Robot
Start MATLAB and run:

text
startup
main_demo
This will:

Initialize the project
Load the robot parameters
Run a multi‑step walking simulation
Display an animation of the RABBIT robot
Controllers
The framework includes several control strategies:

PD Controller

Simple joint‑space proportional‑derivative control used for baseline walking experiments.

Feedback Linearization

Model‑based nonlinear control using the robot dynamics.

Hybrid Zero Dynamics (HZD)

A modern locomotion control framework used to generate stable periodic walking gaits.

Trajectory Optimization
The repository implements direct collocation trajectory optimization for gait generation.

Key components:

dynamic constraints
gait periodicity conditions
optimization cost functions
Important files:

text
direct_collocation.m
gait_constraints.m
optimize_gait.m
rabbit_hzd_trajectory_optimization.m
These allow computation of dynamically feasible walking trajectories.

Visualization
Robot motion is visualized using

text
animate_rabbit.m
Simulation results can also be exported as GIF animations.

Applications
This framework can be used for:

hybrid locomotion research
nonlinear control experiments
trajectory optimization studies
robotics coursework and teaching
rapid testing of walking algorithms
References
E. R. Westervelt, J. W. Grizzle, C. Chevallereau,

J. H. Choi, and B. Morris

Feedback Control of Dynamic Bipedal Robot Locomotion

CRC Press, 2007.

Citation
If you use this project in research or coursework, please cite:

text
@software{rabbit_robot_framework,
  author = {Malek, Mohammad},
  title = {RABBIT Bipedal Robot MATLAB Framework},
  year = {2026},
  url = {https://github.com/M-Malek132/RABBIT-Bipedal-Robot}
}
Author
Mohammad Malek

Robotics & Control Research

Email: Malek.mohammad132@gmail.com

Website: https://malekmohammad.com

GitHub: https://github.com/M-Malek132

License
This project is released under the MIT License.