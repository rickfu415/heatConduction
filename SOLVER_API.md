# 1D Heat Conduction Solver — API Reference for TPS1D-Optimizer

This document describes the solver's public interface as consumed by the web app backend.

---

## Overview

A transient 1D heat conduction solver with adjoint-based thickness optimizer for layered thermal protection systems (TPS).

**Physics:**
- Finite difference, centered in space (2nd order), implicit backward Euler in time (1st order)
- Newton iteration at each timestep (handles nonlinear re-radiation BC)
- Spatially varying material properties across layers
- Re-radiation boundary condition: `q_net = q_in - eps * sigma * (T^4 - T_amb^4)` at `x=0`

**Optimization:**
- Adjoint method: computes `dJ/dk`, `dJ/drho`, `dJ/dcp` per node in one backward pass
- Chain rule maps per-node gradients to layer thickness gradients
- Gradient descent on layer thicknesses; material properties fixed
- Total thickness is free to grow when redistribution alone cannot meet the target

---

## Entry Points

### `optimize_layered(para, max_iter=200, tol=1.0, callback=None)`

**File:** `optimizer.py`

Runs adjoint-based gradient descent to minimize layer thicknesses while keeping backwall temperature below `target`.

#### Input: `para` (pandas.Series)

| Key | Type | Description |
|-----|------|-------------|
| `material function` | str | Must be `'layered'` |
| `numberOfLayers` | int | Number of layers |
| `layerConductivities` | np.array | k per layer (W/m·K) — fixed during optimization |
| `layerDensities` | np.array | rho per layer (kg/m³) — fixed |
| `layerHeatCapacities` | np.array | cp per layer (J/kg·K) — fixed |
| `layerThicknesses` | np.array | Initial layer thicknesses (m) — design variables |
| `length` | float | Sum of `layerThicknesses` (m) |
| `nodesPerLayer` | int | Grid nodes per layer (11 recommended; adjacent layers share 1 boundary node) |
| `duration` | float | Total simulation time (s) |
| `CFL` | float | CFL number for adaptive timestep (default 3.0; increase if Newton doesn't converge) |
| `maxIteration` | int | Max Newton iterations per timestep (50 typical) |
| `convergence` | float | Newton residual convergence threshold (1e-7 typical) |
| `relaxation` | float | Newton underrelaxation factor (0.9; range 0–1, **very sensitive**) |
| `ambientTemperature` | float | Ambient temperature for re-radiation (K) |
| `stefanBoltzmann` | float | Stefan–Boltzmann constant: `5.670374419e-8` W/(m²·K⁴) |
| `emissivity` | float | Surface emissivity at `x=0` |
| `IC value` | float | Uniform initial temperature (K) |
| `x=0 type` | str | `'heatFlux'` or `'fixedTemperature'` |
| `x=0 value` | float or np.array shape (N,2) | Constant flux (W/m²) or time-varying profile `[[t0,q0],[t1,q1],...]` |
| `x=L type` | str | `'heatFlux'` or `'fixedTemperature'` |
| `x=L value` | float | Typically `0.0` (insulated back wall) |
| `re-radiate` | bool | Enable re-radiation at `x=0` |
| `back_wall_temperature_target` | float | Target max backwall temperature (K) |
| `learning_rate` | float | Gradient descent step size (1e-9 typical; thickness gradients are O(1e5)) |
| `print_frequency` | int | Print solver output every N timesteps (default 1) |

#### Output

Returns `(para, history)`:

- **`para`** (pandas.Series): updated with optimized `layerThicknesses` and `length`
- **`history`** (list of dicts): one entry per iteration

| Key | Type | Description |
|-----|------|-------------|
| `iter` | int | Iteration number |
| `T_L` | float | Max backwall temperature this iteration (K) |
| `loss` | float | `0.5 * (T_L - target)^2` |
| `grad_norm` | float | L2 norm of thickness gradient |
| `t_layers` | np.array | Layer thicknesses this iteration (m) |

#### Optional: `callback`

If provided, called each iteration with a dict:

```python
{
    'iter': int,
    'T_L': float,
    'loss': float,
    'grad_norm': float,
    't_layers': list,         # meters
    'total_thickness': float, # meters
}
```

Use this for real-time streaming of progress to the frontend.

---

### `hc.solve(para, verbose=True)`

**File:** `heatConduction.py`

Runs the forward heat conduction simulation.

- Accepts constant or time-varying BCs (2D array `[[t, q], ...]`)
- Returns `(TProfile, cache)`
  - `TProfile`: np.array shape `(numberOfNode, numberOfTimeStep+1)` — temperature field in space and time
  - `cache`: dict containing `TProfile`, `Jacobian`, `flux_history_x0`, `ce_arr`, `cw_arr`, `Log`

---

### `diff.main(para, cache, verbose=False)`

**File:** `differential.py`

Runs the adjoint backward pass on a completed forward solve.

- Loss: `J = 0.5 * (max_t T_backwall(t) - target)^2`
- Returns dict: `{'grad_k', 'grad_rho', 'grad_cp', 'lambda_profile', 'loss'}`

---

## Minimal Usage Example

```python
import numpy as np
import pandas as pd
import parameter
import optimizer

# Build para using parameter.main() or construct manually:
para = parameter.main()

# Override for your problem:
para['layerThicknesses'] = np.array([0.001, 0.002, 0.0005])
para['length'] = float(para['layerThicknesses'].sum())
para['x=0 value'] = np.array([[0, 0], [70, 2e6], [300, 0]])  # time-varying flux
para['back_wall_temperature_target'] = 400.0
para['learning_rate'] = 1e-9

# Run optimizer with streaming callback
def on_progress(d):
    print(f"Iter {d['iter']}: T_L={d['T_L']:.1f}K  total={d['total_thickness']*1000:.2f}mm")

para, history = optimizer.optimize_layered(para, max_iter=100, tol=1.0, callback=on_progress)

print("Optimized thicknesses (mm):", para['layerThicknesses'] * 1000)
print("Total thickness (mm):", para['length'] * 1000)
```

---

## Grid Layout

With `nodesPerLayer = npl` and `n_layers` layers:
- Total nodes: `N = npl * n_layers - (n_layers - 1)`
- Adjacent layers share one boundary node
- `dx` within each layer = `t_layer / (npl - 1)`
- `dx` is non-uniform across layers when thicknesses differ

**Example:** 3 layers, `npl=11` → 31 total nodes (10 intervals per layer)

---

## Adaptive Timestep

`deltaTime` and `numberOfTimeStep` are recomputed from `duration` and `CFL` inside `normalize_conductivity()` at the start of every `hc.solve()` call. The web app can set `CFL` to tune stability vs. speed.

- Per layer: `dt_layer = CFL * dx_layer² / (2 * alpha_layer)`
- `dt = max(dt_candidates)` across layers (limited by slowest layer)
- Capped so there are at least 200 timesteps for transient resolution

---

## Thickness Growth Logic

Before the first gradient step, the optimizer runs a single forward solve to check whether redistribution alone can meet the target — without adding any material.

**Pre-check:** All layers are set to `t_min` except the best insulator (lowest `k`), which receives all remaining thickness. If this best-case redistribution still exceeds the target backwall temperature, growth is enabled.

```
best_insulator = argmin(k_layers)
t_check[best_insulator] = total_thickness - t_min * (n_layers - 1)
T_best = forward_solve(t_check).max_backwall_T

needs_growth = (T_best > target)
```

**During optimization:**

- If `needs_growth = True`: a growth component is added to the gradient each iteration, biased toward thicker low-`k` layers (direction `1/k`, normalized). This breaks the zero-sum symmetry of the chain-rule gradient so total thickness can increase.
- Once `T_L` drops to or below `target`, `needs_growth` is set to `False` permanently — the optimizer switches to redistribution-only mode and total thickness is held constant for the remaining iterations.

**Why this matters:** The chain-rule gradient from the adjoint always sums to zero (moving thickness between layers conserves total). A pure redistribution step can never grow total thickness. The growth component is the only mechanism by which the optimizer can make the TPS stack thicker.

---

## Convergence Tips

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Newton doesn't converge | `relaxation` too high or `CFL` too large | Reduce `relaxation` to 0.7–0.85; reduce `CFL` |
| Optimizer oscillates | `learning_rate` too large | Reduce by 10× |
| Optimizer stalls, T_L barely moves | `learning_rate` too small | Increase by 10× |
| T_L not reaching target after many iters | Total thickness insufficient | Solver will grow thickness automatically; increase `max_iter` |
| Backwall T overshoots target | Loss tolerance `tol` too loose | Reduce `tol` (e.g., 0.1 instead of 1.0) |
