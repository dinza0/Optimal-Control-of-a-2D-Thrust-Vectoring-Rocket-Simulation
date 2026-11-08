# README — `projectpart5_1.m`

## Purpose

`projectpart5_1.m` performs the controllability analysis for the hover and weak time-varying descent linearizations of the variable-mass thrust-vectoring rocket. It computes controllability matrices, singular values, finite-horizon controllability Gramians, Gramian eigenvalues, and a reduced horizontal-subsystem controllability check.

This script corresponds mainly to the early Part I / Section 5.2 controllability work.

## Theory used

The script studies linear systems of the form

```math
\dot{x}=A(t)x+B(t)u.
```

For the hover model, `A_hover` and `B_hover` are constant. The script forms the standard controllability matrix

```math
\mathcal{C}=[B,AB,A^2B,\ldots,A^{n-1}B]
```

and checks its rank and singular values.

For the weak time-varying descent model, the script computes a finite-horizon controllability Gramian. For an LTI or LTV system, the Gramian measures how much control authority exists over a finite time horizon. The script integrates the matrix differential equation

```math
\dot{W}=AW+WA^T+BB^T
```

using `ode45`. Larger Gramian eigenvalues indicate easier-to-control directions; very small eigenvalues indicate weakly controllable directions.

The horizontal subsystem is also isolated because horizontal motion is controlled indirectly through

```text
torque -> angular acceleration -> attitude -> horizontal acceleration -> horizontal position
```

which explains why horizontal controllability is weaker.

## MATLAB requirements

- MATLAB
- No special toolbox required for the main computations beyond standard ODE functions
- Uses `ode45`, `svd`, `eig`, `rank`, and plotting functions

## How to run

1. Open MATLAB.
2. Place `projectpart5_1.m` in your working folder.
3. Run:

```matlab
projectpart5_1
```

The script prints controllability ranks, singular values, condition numbers, Gramian eigenvalues, and creates a semilog plot comparing hover and weak LTV Gramian eigenvalues.

## Main values to change

At the top of the script:

```matlab
p.g     = 9.81;
p.I     = 10;
p.m0    = 20;
p.alpha = 0.005;
tf = 5.0;
n  = 7;
```

Change these to test different physical parameters or time horizons.

### Important tuning variables

- `tf`: finite horizon for Gramian integration. Increasing this usually improves controllability because the system has more time to respond.
- `p.m0`: initial mass. Affects vertical thrust effectiveness.
- `p.I`: moment of inertia. Larger inertia weakens torque authority.
- `p.alpha`: fuel consumption coefficient. Affects the weak descent mass trajectory.

## What the output means

### `rank(C_hover)`
If this equals 7, the hover linearization is controllable in the standard LTI sense.

### `sigma_min(C_hover)`
The smallest singular value of the controllability matrix. Smaller values indicate weak controllability even if rank is full.

### `cond(C_hover)`
The controllability matrix condition number. Large values mean some directions require much more input effort than others.

### Gramian eigenvalues
The finite-horizon Gramian eigenvalues show which state-space directions are easy or hard to influence over the selected horizon. Very small eigenvalues imply weakly controllable directions.

## Expected result

You should see both hover and weak LTV descent report full rank, but with very small minimum Gramian eigenvalues. This means the systems are theoretically controllable, but some directions are difficult to actuate. The weak directions are especially associated with horizontal motion and attitude coupling.

## Troubleshooting

### Rank is not 7
Check that the state ordering and `makeA`, `makeB` functions have not been changed incorrectly.

### Gramian has negative eigenvalues
Small negative values can occur from numerical roundoff. The script symmetrizes the Gramian using

```matlab
W = 0.5*(W + W.');
```

If large negative eigenvalues appear, reduce integration tolerances or check the matrix equations.

### Plot looks flat
The Gramian eigenvalues span many orders of magnitude, so use `semilogy`, as the script already does.
