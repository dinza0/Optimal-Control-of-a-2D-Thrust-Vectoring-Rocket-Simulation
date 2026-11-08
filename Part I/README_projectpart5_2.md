# README — `projectpart5_2.m`

## Purpose

`projectpart5_2.m` solves the finite-horizon LQR regulation problem for both the hover linearization and the weak LTV descent linearization. It corresponds to Part I / Sections 5.3--5.5 of the project.

The script:

1. Defines hover and weak-descent linearized models.
2. Solves the finite-horizon Riccati differential equation backward in time.
3. Simulates the optimal closed-loop system forward in time.
4. Computes state, control, terminal, and total costs.
5. Checks physical thrust, torque, altitude, and dry-mass constraints.
6. Validates PMP/Riccati optimality conditions.
7. Runs a weight tuning sweep.
8. Runs initial-condition sensitivity cases.
9. Generates plots for states, controls, Hamiltonian, and trajectories.

## Theory used

The finite-horizon LQR problem is

```math
J=\int_0^{t_f}\left(x^TQx+u^TRu\right)dt+x(t_f)^TQ_fx(t_f)
```

subject to

```math
\dot{x}=A(t)x+B(t)u.
```

The Hamiltonian is

```math
H=x^TQx+u^TRu+\lambda^T(Ax+Bu).
```

PMP gives

```math
\dot{x}=Ax+Bu,
```

```math
\dot{\lambda}=-2Qx-A^T\lambda,
```

```math
2Ru+B^T\lambda=0.
```

For an LQR problem, the costate is written as

```math
\lambda(t)=2P(t)x(t),
```

which gives the feedback law

```math
u^*(t)=-R^{-1}B(t)^TP(t)x(t).
```

The Riccati equation solved by the script is

```math
-\dot{P}=A^TP+PA-PBR^{-1}B^TP+Q,
```

with terminal condition

```math
P(t_f)=Q_f.
```

Because the dynamics are linear and the cost is quadratic with `R > 0`, satisfying the Riccati/PMP conditions gives the unconstrained global optimum. Physical actuator feasibility is then checked after the solution is computed.

## MATLAB requirements

- MATLAB
- Uses `ode45`
- Uses standard plotting functions
- No optimization toolbox is needed for the core Riccati simulation

## How to run

1. Open MATLAB.
2. Place `projectpart5_2.m` in your working folder.
3. Run:

```matlab
projectpart5_2
```

The script automatically solves both the hover and weak LTV descent cases, prints numerical summaries, runs tuning/sensitivity analyses, and generates plots.

## Main configuration section

Most values are changed in the first section labeled `User Configuration`.

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
```

## Changing the initial condition

The initial perturbation state is

```matlab
x0 = [8; -4; 2; -1; deg2rad(5); deg2rad(2); 0];
```

The state order is

```text
[dpx; dpz; dvx; dvz; dtheta; domega; dm]
```

Change these entries to test different regulation scenarios.

## Changing LQR weights

The tuned weights are

```matlab
Q  = diag([10, 10, 1, 1, 50, 5, 0.1]);
R  = diag([0.05, 2.7]);
Qf = diag([100, 100, 10, 10, 100, 20, 1]);
```

### What each matrix does

- `Q`: running state error penalty
- `R`: running control effort penalty
- `Qf`: terminal state error penalty

### Common tuning changes

- Increase `R(2,2)` to reduce torque usage.
- Increase `R(1,1)` to reduce thrust perturbation usage.
- Increase attitude weight `Q(5,5)` to penalize pitch angle more.
- Increase terminal weights in `Qf` to reduce final state error.
- Decrease `R` if the controller is too weak, but this can violate actuator limits.

## Weight tuning sweep

The section labeled `5.5 Weight Tuning Sweep` tests multiple sets of `Q`, `R`, and `Qf`. Add another case using the same pattern:

```matlab
weightCases{end+1}.name = 'My test case';
weightCases{end}.Q  = diag([...]);
weightCases{end}.R  = diag([...]);
weightCases{end}.Qf = diag([...]);
```

## Sensitivity cases

The section labeled `Initial Condition Sensitivity Study` changes the initial state to test robustness. Add a case by appending to `icCases`.

## Solver details

### Riccati equation
Solved backward using:

```matlab
ode45(riccati_rhs, [cfg.tf 0], p0_vec, opts)
```

### Closed-loop dynamics
Simulated forward using:

```matlab
ode45(closed_loop_rhs, [0 cfg.tf], x0, opts)
```

## Output interpretation

### `J_state`
The accumulated state error cost.

### `J_control`
The accumulated control effort cost.

### `J_final`
The terminal penalty.

### `J_total`
Total objective value.

### `max |delta T|`, `max |tau|`
Indicate how aggressively the controller uses inputs.

### Physical thrust bounds
The script converts perturbation thrust to physical thrust using

```matlab
Tphys = cfg.m_hover*cfg.g + dT;
```

and checks whether it lies inside `[Tmin, Tmax]`.

### PMP validation
The stationarity residual

```math
\|2Ru+B^T\lambda\|
```

should be near zero. The terminal costate condition

```math
\|\lambda(t_f)-2Q_fx(t_f)\|
```

should also be near zero.

## Expected result

The final tuned weights should satisfy thrust and torque limits for both hover and weak LTV descent. The terminal state norm will not be zero because control effort is penalized and actuator feasibility required less aggressive tuning.

## Troubleshooting

### Torque bound violated
Increase `R(2,2)`.

### Thrust bound violated
Increase `R(1,1)` or reduce vertical/terminal weights.

### Final error too large
Increase `Qf` or relevant entries of `Q`, but check actuator limits again.

### Plots look identical for hover and weak descent
This can happen for short horizons and small mass variation. Increase `tf` or `alpha` to highlight LTV effects.
