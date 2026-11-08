# README — `projectpart7b.m`

## Purpose

`projectpart7b.m` solves Part III Phase B: return-and-land energy-optimal control. It computes two landing solutions:

1. Modified ground landing manifold `M1x`
2. Circular platform landing manifold `M2`

The script uses PMP boundary conditions and MATLAB's `bvp4c` solver.

## Theory used

Phase B minimizes input energy:

```math
J=\int u^TRu\,dt.
```

A small time regularization is included numerically:

```math
J=\int (u^TRu+\rho)\,dt.
```

The Hamiltonian is

```math
H=u^TRu+\rho+\lambda^T(Ax+Bu).
```

PMP gives

```math
\dot{x}=Ax+Bu,
```

```math
\dot{\lambda}=-A^T\lambda,
```

```math
2Ru+B^T\lambda=0,
```

so

```math
u^*=-\frac{1}{2}R^{-1}B^T\lambda.
```

Because final time is free,

```math
H(t_f)=0.
```

Because the endpoint lies on a terminal manifold,

```math
\lambda(t_f)=D\psi(x_f)^T\mu.
```

## MATLAB requirements

- MATLAB
- Uses `bvp4c`
- Uses standard plotting functions

## How to run

1. Open MATLAB.
2. Put `projectpart7b.m` in the working folder.
3. Run:

```matlab
projectpart7b
```

The script solves `M1x`, solves `M2`, prints a comparison table, and plots landing trajectories, states, controls, Hamiltonians, thrust vectors, and the circular platform graphic.

## Important note about M1

The original project manifold was

```math
M_1=\{p_z=0,\ v_z=0,\ \theta=0\}.
```

This leaves final horizontal position free. With pure energy minimization and free final time, this caused numerical degeneracy and trivial slow-descent solutions.

The code uses a modified point-landing version:

```math
M_{1x}=\{p_x=p_{x,land},\ p_z=0,\ v_z=0,\ \theta=0\}.
```

This should be stated clearly in the report.

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

### Energy weights

```matlab
p.R = diag([0.05, 2.7]);
p.rho = 1e-2;
```

- Increase `p.R(1,1)` to reduce thrust perturbation.
- Increase `p.R(2,2)` to reduce torque.
- Increase `p.rho` to penalize long final times more strongly.

### Circular platform

```matlab
p.pc = 0;
p.r  = 2;
```

- `p.pc`: horizontal center of the platform
- `p.r`: platform radius

### Phase B initial sky state

```matlab
x_sky = [10; 25; 0; 0; 0; 0; dm_sky];
```

State order:

```text
[px; pz; vx; vz; theta; omega; dm]
```

Change this to match the final state of Phase A if desired.

### Ground landing x-location

```matlab
p.land_px = 0;
```

This is the desired final horizontal position for the modified ground landing case.

## Solver details

The script uses a normalized time variable:

```math
s\in[0,1],\quad t=s t_f.
```

To ensure positive final time, it parameterizes:

```math
t_f=e^\eta.
```

The unknown BVP parameters are:

- `eta`, which determines final time
- terminal manifold multipliers `mu`

The boundary conditions enforce:

1. Initial state equals `x_sky`
2. Terminal manifold constraints
3. Terminal costate transversality
4. Free-final-time Hamiltonian condition

## Landing manifolds in code

### `M1x`

```matlab
psi = [xf(1)-p.land_px; xf(2); xf(4); xf(5)];
```

This enforces:

```text
px = p.land_px
pz = 0
vz = 0
theta = 0
```

### `M2`

```matlab
psi = [(xf(1)-p.pc)^2 + xf(2)^2 - p.r^2;
       xf(4);
       xf(5)];
```

This enforces:

```text
final position on circle
vz = 0
theta = 0
```

## Output interpretation

### `tf`
Final landing time.

### `J`
Energy/time-regularized cost.

### `H(tf)`
Should be near zero for a free-final-time problem.

### `Hamiltonian variation`
Should be small if the BVP solution is consistent.

### `terminal manifold residual`
Should be near zero.

### Controls
For an energy-optimal problem, controls should generally be smooth and low magnitude, unlike Phase A bang-bang controls.

## Expected result

Both `M1x` and `M2` should solve with positive final time, low thrust variations around hover thrust, very small torque, and near-zero terminal manifold residuals.

## Troubleshooting

### Negative final time
This should be prevented by `tf = exp(eta)`. Make sure that transformation is still in the code.

### `bvp4c` singular Jacobian
Try one or more of the following:

- Increase `p.rho`
- Improve `tf_guess`
- Improve `xf_guess`
- Use `M1x` instead of original `M1`
- Reduce the difficulty of `x_sky`

### Circular platform does not appear correctly
Check:

```matlab
p.pc
p.r
```

and the platform plotting section.

### Torque is almost zero
This is not necessarily an error. For an energy-minimizing landing with small attitude requirements, the optimal solution may avoid torque almost entirely.
