# README — `projectpart7a.m`

## Purpose

`projectpart7a.m` solves Part III Phase A: minimum-time ascent to a sky rendezvous state. It first solves a smoothed indirect shooting problem using an epsilon continuation method, then refines the result into an exact bang-bang switch-time parameterization.

This script is used to obtain the final Phase A exact bang-bang trajectory, thrust/torque plots, Hamiltonian plot, and switching function plots.

## Theory used

Phase A is a minimum-time problem:

```math
J=t_f=\int_0^{t_f}1\,dt.
```

The Hamiltonian is

```math
H=1+\lambda^T(Ax+Bu).
```

Because the controls enter the dynamics linearly and are bounded, the Hamiltonian is affine in the controls. PMP therefore implies bang-bang control unless a singular arc occurs.

The switching functions are the coefficients of the controls in the Hamiltonian:

```math
\phi_T = \frac{1}{m_h}\lambda_{v_z}-\alpha\lambda_m,
```

```math
\phi_\tau = \frac{1}{I}\lambda_\omega.
```

The signs of these switching functions determine whether the optimal control chooses its upper or lower bound.

The free-final-time condition is

```math
H(t_f)=0.
```

The mass is free at the terminal state, so

```math
\lambda_m(t_f)=0.
```

## MATLAB requirements

- MATLAB
- Optimization Toolbox required for `lsqnonlin`
- Uses `ode45`

## How to run

1. Open MATLAB.
2. Put `projectpart7a.m` in the working folder.
3. Run:

```matlab
projectpart7a
```

The script will:

1. Solve a smoothed minimum-time shooting problem for decreasing epsilon values.
2. Choose the best smoothed solution.
3. Simulate the smoothed solution.
4. Reconstruct approximate bang-bang controls from switching functions.
5. Solve an exact bang-bang switch-time problem.
6. Plot exact bang-bang trajectory, controls, Hamiltonian, and switching functions.

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

### Initial state

```matlab
x0 = [0;0;0;0;0;0;0];
```

State order:

```text
[px; pz; vx; vz; theta; omega; dm]
```

### Terminal target

The active target is

```matlab
xsky6 = [10;15;0;0;0;0];
```

Only the first six states are constrained. Mass is free.

To test another target, change `xsky6`, for example:

```matlab
xsky6 = [20;25;0;0;0;0];
xsky6 = [0;10;0;0;0;0];
xsky6 = [10;15;0;2;0;0];
```

## Epsilon continuation

The smoothed solve uses:

```matlab
eps_list = [20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002];
```

Smaller epsilon values make the control closer to bang-bang but also make the shooting problem harder and more singular.

If convergence fails, remove the smallest epsilon values or improve the initial guess.

## Initial guess

The initial shooting guess is

```matlab
zguess = [-1; -1; -0.5; -1; 0.1; 0.1; 0.1; 6.0];
```

This vector contains:

```text
[lambda0(1:7); tf]
```

Change the last entry to adjust the initial final-time guess.

## Solver details

### Smoothed solve
Uses `lsqnonlin` with bounds on costates and final time:

```matlab
lb = [-1e3*ones(7,1); 0.5];
ub = [ 1e3*ones(7,1); 30];
```

The residual contains terminal state error, mass-free transversality, and Hamiltonian terminal condition.

### Exact bang-bang switch-time solve
After the smoothed solution identifies the switching structure, the script solves for:

```text
[lambda0(1:7); tf; thrust switch times; torque switch times]
```

using `lsqnonlin`.

The final implemented structure is:

```text
T:   Tmax -> Tmin -> Tmax
tau: taumax -> taumin -> taumax -> taumin
```

## Changing switch structure

Look for the section labeled:

```matlab
Phase A: Exact Bang-Bang Switch-Time Parameterization
```

The initial guesses are:

```matlab
tT10 = 1.90;
tT20 = 2.80;
ttau10 = 0.48;
ttau20 = 1.62;
ttau30 = 2.68;
```

If you change the terminal target, the switching structure may change. Use the smoothed solution plots to estimate new switch times.

## Output interpretation

### `terminal residual norm`
Should be very small, ideally less than `1e-6`. Your final exact solution reached about `1e-13`.

### `Hamiltonian final`
Should be near zero for a free-final-time problem.

### `Hamiltonian variation`
Should be very small for the exact bang-bang solution.

### `max physical thrust`, `min physical thrust`
Should hit `300` and `0` for the exact bang-bang result.

### `max abs torque`
Should hit `30` if torque saturates.

## Which plots to use in the report

Use the exact bang-bang plots as final results:

- Exact bang-bang trajectory
- Exact physical thrust and torque
- Exact Hamiltonian
- Switching functions

The smoothed plots are useful for methodology but are not the final solution.

## Troubleshooting

### `lsqnonlin` fails for small epsilon
This is expected. The problem becomes nearly discontinuous. Use the smallest epsilon that converges reliably, then switch to switch-time parameterization.

### Final residual is large in exact bang-bang solve
The assumed switch structure may be too restrictive. Add or remove switch times based on the smoothed switching functions.

### Hamiltonian is not near zero
Check the residual scaling, final-time condition, and whether the switch times are ordered correctly.

### Exact solve finds local minimum possible
This can still be acceptable if terminal residual, Hamiltonian residual, and switching functions are all very small.
