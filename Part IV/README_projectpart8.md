# README — `projectpart8.m`

## Purpose

`projectpart8.m` solves Part IV: nonlinear optimal control near hover. It compares three cases:

1. Linear optimal controller on linear dynamics
2. Same linear controller on nonlinear dynamics
3. Nonlinear optimal controller on nonlinear dynamics

The script is used to determine how accurate the linearized controller remains near hover and how much the full nonlinear optimal solution improves performance.

## Theory used

The full nonlinear variable-mass rocket dynamics are

```math
\dot{p}_x=v_x,
```

```math
\dot{p}_z=v_z,
```

```math
\dot{v}_x=\frac{T}{m}\sin\theta,
```

```math
\dot{v}_z=\frac{T}{m}\cos\theta-g,
```

```math
\dot{\theta}=\omega,
```

```math
\dot{\omega}=\frac{\tau}{I},
```

```math
\dot{m}=-\alpha T.
```

The cost is

```math
J=\int_0^{t_f}(e^TQe+u^TRu)dt,
```

where

```math
e=x-x_{hover}.
```

The nonlinear Hamiltonian is

```math
H=e^TQe+u^TRu+\lambda^Tf(x,u).
```

The nonlinear optimal control satisfies stationarity:

```math
\frac{\partial H}{\partial u}=0.
```

The terminal state is free and there is no terminal cost, so

```math
\lambda(t_f)=0.
```

Since final time is fixed, the condition `H(tf)=0` is not required.

## MATLAB requirements

- MATLAB
- Optimization Toolbox required for `lsqnonlin`
- Uses `ode45`

## How to run

1. Open MATLAB.
2. Place `projectpart8.m` in your working folder.
3. Run:

```matlab
projectpart8
```

The script will solve all three comparison cases, print a summary table, and create comparison plots for states, controls, trajectory, and Hamiltonian.

## Main values to change

### Physical parameters

```matlab
p.g = 9.81;
p.I = 10;
p.mh = 20;
p.alpha = 0.005;
p.Tmin = 0;
p.Tmax = 300;
p.taumin = -30;
p.taumax = 30;
```

### Final time

```matlab
p.tf = 2.5;
```

Part IV is intended for short near-hover maneuvers. Increasing `p.tf` can make nonlinear shooting harder.

### Hover state

```matlab
p.xhover = [0; 50; 0; 0; 0; 0; p.mh];
```

This is the state about which the linear controller is designed.

### Weights

```matlab
p.Q = diag([10, 10, 1, 1, 50, 5, 0.1]);
p.R = diag([0.05, 2.7]);
p.Qf = zeros(7);
```

No terminal cost is used in Part IV.

### Initial condition

```matlab
e0 = [1.5; -1; 0.3; -0.2; deg2rad(2); deg2rad(0.5); 0];
x0 = p.xhover + e0;
```

State order:

```text
[px; pz; vx; vz; theta; omega; m]
```

`e0` is the perturbation from hover.

## Solver details

### Case 1: Linear optimal on linear dynamics
Uses finite-horizon Riccati equation solved backward with `ode45`, then simulates the closed-loop linear system forward.

### Case 2: Linear controller on nonlinear dynamics
Uses the same Riccati feedback controller from Case 1, but applies it to the full nonlinear dynamics.

### Case 3: Nonlinear optimal on nonlinear dynamics
Uses indirect shooting with `lsqnonlin`.

Unknown variable:

```text
lambda0 = initial costate
```

Residual:

```text
lambda(tf)
```

The initial guess is generated from the linear Riccati solution:

```matlab
lambda0_guess = 2*P0*e0;
```

This is important for convergence.

## Output interpretation

### `J`
Total running cost.

### `||e(tf)||`
Final state error norm. This is not expected to be zero because there is no terminal constraint or terminal cost.

### `maxT`
Maximum physical thrust used.

### `maxTau`
Maximum torque magnitude.

### `Hvar`
Hamiltonian variation.

For fixed-final-time autonomous optimal solutions, the Hamiltonian should be approximately constant but does not have to be zero.

### Nonlinear terminal costate residual

```math
\|\lambda(t_f)\|
```

Should be near zero for a valid nonlinear optimal solution.

### Nonlinear stationarity residual

```math
\max_t\|H_u\|
```

Should be near zero.

## Expected result

The linear and nonlinear trajectories should be close for the near-hover initial condition. The nonlinear optimal solution should have slightly lower cost than applying the linear controller to the nonlinear dynamics.

## Important interpretation

- Case 1 is optimal for the linearized model.
- Case 2 is not nonlinear-optimal, but tests whether the linear controller remains viable on the nonlinear plant.
- Case 3 is the nonlinear PMP solution and should be judged by terminal costate residual, stationarity residual, and Hamiltonian consistency.

## Troubleshooting

### Nonlinear shooting fails or residual is large
Try:

- Reduce `e0`
- Reduce `p.tf`
- Increase `MaxIterations` and `MaxFunctionEvaluations`
- Use a better `lambda0_guess`
- Avoid large initial attitude angles

### Nonlinear optimal trajectory looks unreasonable
The initial condition may be outside the near-hover validity region, or the shooting method may have converged to a poor local solution.

### Hamiltonian is not zero
That is okay for Part IV because final time is fixed. Check Hamiltonian variation instead.

### Final error is not zero
That is expected because terminal state is free and `Qf = zeros(7)`.

### Controls saturate
For Part IV near-hover, saturation usually indicates the initial condition is too large or weights are too aggressive.
