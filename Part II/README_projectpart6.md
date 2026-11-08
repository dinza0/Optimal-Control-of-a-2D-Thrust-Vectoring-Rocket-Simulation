# README — `projectpart6.m`

## Purpose

`projectpart6.m` solves Part II of the project: optimal trajectory tracking. It first generates a Part I optimal reference trajectory, then designs an affine finite-horizon tracking controller that drives a different initial condition toward that reference trajectory.

The script also performs a sensitivity analysis by testing larger and altered initial conditions.

## Theory used

The tracking problem is

```math
J=\int_0^{t_f}\left[(x-x_{ref})^TQ(x-x_{ref})+u^TRu\right]dt
```

subject to

```math
\dot{x}=A(t)x+B(t)u.
```

The Hamiltonian is

```math
H=(x-x_{ref})^TQ(x-x_{ref})+u^TRu+\lambda^T(Ax+Bu).
```

PMP gives

```math
\dot{\lambda}=-2Q(x-x_{ref})-A^T\lambda,
```

and stationarity gives

```math
2Ru+B^T\lambda=0.
```

For tracking, the costate is written as

```math
\lambda(t)=2(P(t)x(t)+s(t)).
```

The resulting affine control law is

```math
u^*(t)=-R^{-1}B(t)^T(P(t)x(t)+s(t)).
```

The matrix `P(t)` satisfies the usual Riccati equation, while the vector `s(t)` accounts for the nonzero reference trajectory. Since the Part II project statement does not include a terminal tracking cost, the terminal conditions are

```math
P(t_f)=0,
```

and

```math
s(t_f)=0.
```

This also means the final tracking error is not forced to zero.

## MATLAB requirements

- MATLAB
- Uses `ode45`
- Uses standard plotting functions

## How to run

1. Open MATLAB.
2. Put `projectpart6.m` in the working folder.
3. Run:

```matlab
projectpart6
```

The script will:

1. Generate a Part I weak LTV descent reference trajectory.
2. Choose a different tracking initial condition.
3. Solve the affine tracking Riccati equations.
4. Simulate the tracking system.
5. Print a tracking summary.
6. Run sensitivity cases.
7. Generate state/reference, error, control, Hamiltonian, and trajectory plots.

## Main configuration values

At the top of the script:

```matlab
cfg.g      = 9.81;
cfg.I      = 10;
cfg.m0     = 20;
cfg.alpha  = 0.005;
cfg.Tmin   = 0;
cfg.Tmax   = 300;
cfg.taumin = -30;
cfg.taumax = 30;
cfg.mdry   = 5;
cfg.tf     = 5.0;
cfg.Nplot  = 1001;
```

## Changing the reference trajectory

The Part I reference initial condition is

```matlab
x0_partI = [8; -4; 2; -1; deg2rad(5); deg2rad(2); 0];
```

This is used to generate `xref_grid`. Change `x0_partI` if you want a different reference trajectory.

## Changing the tracking initial condition

The tracking initial condition is

```matlab
x0_track = [-6; 6; -3; 2; deg2rad(-8); deg2rad(-3); 0];
```

State order:

```text
[dpx; dpz; dvx; dvz; dtheta; domega; dm]
```

Increase or decrease these values to test harder or easier tracking problems.

## Changing weights

The current tracking weights are

```matlab
Q  = diag([10, 10, 1, 1, 50, 5, 0.1]);
R  = diag([0.04, 7.5]);
Qf = diag([100, 100, 10, 10, 100, 20, 1]);
```

Note: `Qf` is used in generating the Part I reference trajectory, but the actual tracking problem does not include a terminal tracking cost.

### Tuning recommendations

- Increase `R(2,2)` if torque violates bounds.
- Increase `R(1,1)` if thrust violates bounds.
- Increase entries of `Q` to improve tracking of specific states.
- To reduce final tracking error, a terminal cost would be needed, but that is not part of the assigned Part II formulation.

## Sensitivity analysis

The script tests these cases:

- Nominal tracking initial condition
- 25 percent larger initial condition
- 50 percent larger initial condition
- High horizontal velocity error
- High attitude error

To add another case, append to `icCases`:

```matlab
icCases{end+1}.name = 'My case';
icCases{end}.x0 = [...];
```

## Solver details

### Part I reference solve
Uses finite-horizon Riccati integration backward with `ode45`, then closed-loop simulation forward with `ode45`.

### Tracking solve
The tracking Riccati variables `P(t)` and `s(t)` are stacked into a single vector and integrated backward using `ode45`.

### Tracking closed-loop simulation
The affine control law is evaluated during forward simulation using `ode45`.

## Output interpretation

### `J_state`
Integral of tracking error cost.

### `J_control`
Integral of control effort cost.

### `J_total`
Total tracking objective.

### `Max tracking error norm`
Largest value of `norm(x - xref)` over the trajectory.

### `Final tracking error norm`
Tracking error at final time. This is not expected to be zero unless a terminal cost or terminal constraint is added.

### PMP validation
The script checks

```math
\max_t\|2Ru+B^T\lambda\|
```

and

```math
\|\lambda(t_f)\|.
```

Both should be near zero for the stated tracking problem.

## Expected result

The nominal case should be feasible, with thrust and torque inside physical bounds. Sensitivity cases may become infeasible because larger deviations require more torque authority.

## Important interpretation

The controller is optimal for the stated running-cost tracking problem, but it does not guarantee zero terminal error. To force arrival at the final reference state, add a terminal cost or terminal equality constraint.

## Troubleshooting

### Torque bound violated
Increase `R(2,2)` or reduce initial attitude/horizontal velocity error.

### Final error is large
This is expected without terminal cost. Add terminal cost only if allowed by the project instructions.

### Sensitivity cases infeasible
This indicates the affine tracking controller is locally valid but not robust to large deviations.
