# Optimal Control of a Variable-Mass Thrust-Vectoring Rocket

This repository contains MATLAB code and report files for an optimal control project involving a simplified planar variable-mass thrust-vectoring rocket. The project develops and compares linearized LQR regulation, trajectory tracking, minimum-time bang-bang control, energy-optimal landing, and nonlinear near-hover optimal control.

## Project Overview

The rocket state is modeled as

\[
x = [p_x,\ p_z,\ v_x,\ v_z,\ \theta,\ \omega,\ m]^T
\]

with thrust and torque controls. The project studies how different optimal control formulations affect the resulting trajectory, control effort, feasibility, and optimality conditions.

Main topics include:

- Hover and weak-descent linearization
- Controllability analysis
- Finite-horizon LQR regulation
- Optimal trajectory tracking
- Pontryagin Minimum Principle
- Bang-bang minimum-time control
- Switching functions
- Free-final-time transversality conditions
- Energy-optimal landing
- Nonlinear optimal control near hover

## Repository Contents

| File | Description |
|---|---|
| `projectpart5_1.m` | Hover linearization, controllability, and Part I setup |
| `projectpart5_2.m` | Weak-descent LTV model and finite-horizon LQR regulation |
| `projectpart6.m` | Optimal trajectory tracking using an affine Riccati tracking law |
| `projectpart7a.m` | Phase A minimum-time ascent and exact bang-bang solution |
| `projectpart7b.m` | Phase B energy-optimal return and landing problem |
| `projectpart8.m` | Nonlinear near-hover optimal control comparison |
| `README_projectpart*.md` | Individual documentation files for each MATLAB script |
| `report/` | Final project report, figures, and supporting files |

## Requirements

The code was developed in MATLAB and uses common MATLAB numerical tools, including:

- `ode45`
- `ode113`
- `fsolve`
- `lsqnonlin`
- `bvp4c`
- MATLAB plotting utilities

Some scripts require the Optimization Toolbox.

## How to Run

Open MATLAB, place all project scripts in the same folder, and run each script directly from the command window or editor.

Example:

```matlab
run projectpart5_1.m
run projectpart5_2.m
run projectpart6.m
run projectpart7a.m
run projectpart7b.m
run projectpart8.m
````

Each script generates plots, numerical summaries, and optimality checks relevant to its project section.

## Changing Parameters

Most important parameters are defined near the top of each script, including:

```matlab
g = 9.81;
mh = 20;
I = 10;
alpha = 0.005;
Tmin = 0;
Tmax = 300;
taumin = -30;
taumax = 30;
```

Controller weights such as `Q`, `R`, and `Qf` can also be changed inside the corresponding script to tune tracking accuracy, control effort, and actuator feasibility.

## Notes

The linearized problems use thrust perturbation,

[
\delta T = T - m_h g
]

while physical feasibility is checked using the reconstructed physical thrust,

[
T = m_h g + \delta T.
]

The Phase A bang-bang solution uses the linearized model with physical actuator bounds. The nonlinear model is studied separately in Part IV.

## Author

Aman D.
UCLA MAE 270C Optimal Control Project

```
```
