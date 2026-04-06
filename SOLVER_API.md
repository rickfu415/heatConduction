# 1D Heat Conduction Solver — API Reference for TPS1D-Optimizer

This document describes the solver's public interface as consumed by the web app backend.

---

## Overview

A transient 1D heat conduction solver with adjoint-based thickness optimizer for layered thermal protection systems (TPS). Optionally coupled to `MaterialStateSolver/material_engine` for evolving material properties (pyrolysis, char formation, T-dependent k/ρ/cp) and volumetric reaction heat source Q.

**Physics:**
- Finite difference, centered in space (2nd order), implicit backward Euler in time (1st order)
- Newton iteration at each timestep (handles nonlinear re-radiation BC)
- Spatially varying material properties across layers
- Re-radiation boundary condition: `q_net = q_in - eps * sigma * (T^4 - T_amb^4)` at `x=0`
- Optional volumetric heat source Q (W/m³) in the PDE: `ρ·cp·∂T/∂t = ∂/∂x(k·∂T/∂x) + Q`
- Optional material state coupling: per-layer `MaterialEngine` evolves k, ρ, cp, Q each timestep based on temperature history (pyrolysis, char, porosity, effective-medium conductivity)

**Optimization:**
- Adjoint method: computes `dJ/dk`, `dJ/drho`, `dJ/dcp` per node in one backward pass (constant-property only)
- Chain rule maps per-node gradients to layer thickness gradients
- Gradient descent on layer thicknesses; material properties fixed
- SLSQP mass minimizer: finite-difference Jacobians, compatible with evolving properties
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

### `hc.solve(para, verbose=True, material_hook=None)`

**File:** `heatConduction.py`

Runs the forward heat conduction simulation.

- Accepts constant or time-varying BCs (2D array `[[t, q], ...]`)
- Optional `material_hook(para, cache, timeStep)` is called after each converged Newton step (before `storeUpdateResult`), allowing per-timestep updates to `para['conductivity'/'density'/'heatCapacity'/'volumetricHeatSource']`. See [Material Coupling](#material-coupling) below.
- Returns `(TProfile, cache)`
  - `TProfile`: np.array shape `(numberOfNode, numberOfTimeStep+1)` — temperature field in space and time
  - `cache`: dict containing `TProfile`, `Jacobian`, `flux_history_x0`, `ce_arr`, `cw_arr`, `Log`

---

### `diff.main(para, cache, verbose=False)`

**File:** `differential.py`

Runs the adjoint backward pass on a completed forward solve.

- Loss: `J = 0.5 * (max_t T_backwall(t) - target)^2`
- Returns dict: `{'grad_k', 'grad_rho', 'grad_cp', 'lambda_profile', 'loss'}`

> **NOTE:** The adjoint assumes time-invariant k, ρ, cp and no volumetric source Q. It is NOT valid when `hc.solve` was called with a `material_hook`. Use `optimize_mass_slsqp` (FD-based) for evolving-property optimization.

---

### `optimizer.optimize_mass_slsqp(para_base, t0, T_bw_limit, layer_service_temps, t_min, ...)`

**File:** `optimizer.py`

Constrained areal-mass minimization via SLSQP. Compatible with evolving material properties.

```
minimize    m(t) = sum_i rho_i * t_i
subject to  T_backwall_max(t)  <= T_bw_limit
            T_layer_i_max(t)   <= layer_service_temps[i]   (if finite)
            t_i >= t_min[i]
```

#### Key Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `para_base` | pd.Series | Base parameter set (never mutated; copies made internally) |
| `t0` | np.array (n_layers,) | Initial layer thicknesses (m) |
| `T_bw_limit` | float | Max allowed backwall temperature (K) |
| `layer_service_temps` | list[float] | Per-layer service T limits (K); `np.inf` to skip |
| `t_min` | float or np.array | Lower bound(s) on thicknesses (m) |
| `t_max` | float, np.array, or None | Upper bound(s); None = no upper bound |
| `max_iter` | int | Max SLSQP iterations (default 300) |
| `callback` | callable or None | Called each iteration with progress dict |

#### Enabling Material Coupling

Set `para_base['material_engine_yamls']` to a list of material YAML paths or names (one per layer) before calling. Each forward solve inside the optimizer builds a fresh `LayeredMaterialCoupler` so FD perturbations start from initial material state.

```python
para['material_engine_yamls'] = ['pica', 'steel_304']
res = optimizer.optimize_mass_slsqp(para_base=para, ...)
```

#### Output: `scipy.optimize.OptimizeResult`

Standard scipy result plus:
- `.history`: list of per-iteration dicts (`mass_kg_m2`, `t_layers_mm`, `T_bw_max`, etc.)
- `.para_opt`: updated para Series at the optimum

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

---

## Material Coupling

**File:** `material_coupling.py`

Optional coupling to `MaterialStateSolver/material_engine` for temperature-dependent, evolving material properties. The coupling is **opt-in** and **lazy-imported** — `material_engine` is only loaded when coupling is actually used. The solver runs standalone without it.

### Architecture

Each TPS layer gets its own `MaterialEngine` (loaded from a YAML kinetic model) and `MaterialState` (per-node component densities, structural variables). A `material_hook` callback is called after each converged Newton step inside `hc.solve()`:

```
for timeStep in 1..N:
    newtonIteration()        # solve heat equation with current k, rho, cp, Q
    material_hook()          # evolve state, write updated k, rho, cp, Q into para
    storeUpdateResult()      # save T to TProfile, advance T0
```

At the hook call, `cache['T0']` = T_old (previous step) and `cache['T']` = T_new (just converged). The hook calls `MaterialEngine.material_update(T_old, T_new, state, dt)` per layer, then writes per-node arrays back into `para['conductivity'/'density'/'heatCapacity'/'volumetricHeatSource']`.

### Key Classes and Functions

#### `resolve_material_yaml(name_or_path) -> Path`

Resolves a material name (e.g., `'pica'`) or path to an absolute YAML path. Search order:
1. Literal path (absolute or cwd-relative).
2. `$MATERIAL_STATE_SOLVER_PATH/tps_material_db/models/**/<name>.yaml`
3. Default: `../../MaterialStateSolver/tps_material_db/models/` (relative to `material_coupling.py`).

#### `LayeredMaterialCoupler(para, yaml_paths, *, zero_Q=False)`

Creates one `MaterialEngine` + `MaterialState` per layer, sized to `nodesPerLayer` nodes each.

| Method | Description |
|--------|-------------|
| `hook(para, cache, timeStep)` | The callable passed to `hc.solve(material_hook=...)` |
| `reset(para)` | Re-create states from initial YAML and re-seed properties |
| `states` | List of `MaterialState` objects (one per layer); inspect for component densities, porosity, etc. |

#### `make_layered_coupler(para, yaml_paths=None, *, zero_Q=False)`

Factory. If `yaml_paths` is None, resolves from `para['materials']`.

### Volumetric Heat Source Q

When coupling is active, the `MaterialEngine` returns a per-node Q (W/m³) from reaction enthalpies (e.g., endothermic pyrolysis for PICA: Q < 0). This is added to the PDE residual in `assemble()`:

```
F = T - T0 - dt/(rho*cp) * diffusion - dt/(rho*cp) * Q
```

Q is treated as frozen within each Newton iteration (`dQ/dT = 0` in the Jacobian). It is refreshed by the hook between timesteps. For stiff reactions, reduce `CFL` to compensate.

Set `zero_Q=True` in the coupler to disable the source term while still evolving k/ρ/cp (useful for debugging).

### Shared Boundary Nodes

Adjacent layers share a boundary node. When writing per-node properties, layers are iterated in order and **last-write-wins** — the shared node gets layer l+1's properties. This matches `parameter.normalize_conductivity()`.

### Usage Examples

#### Forward solve with coupling

```python
import heatConduction as hc
from material_coupling import make_layered_coupler

coupler = make_layered_coupler(para, ['pica', 'steel_304'])
TProfile, cache = hc.solve(para, material_hook=coupler.hook)

# Inspect final material state
print(coupler.states[0].component_densities['phenolic_resin'])  # PICA layer
print(coupler.states[0].structural_variables['porosity'])
```

#### SLSQP optimizer with coupling

```python
import optimizer

para['material_engine_yamls'] = ['pica', 'steel_304']
res = optimizer.optimize_mass_slsqp(
    para_base=para, t0=t0,
    T_bw_limit=450.0,
    layer_service_temps=[np.inf, np.inf],
    t_min=np.full(2, 0.001),
)
```

#### Standalone script

```python
python run_coupled_all_cases.py   # sys.path hack included, no pip install needed
```

### Adjoint Compatibility

| Optimizer | Compatible with coupling? | Notes |
|-----------|--------------------------|-------|
| `optimize_mass_slsqp` | Yes | FD Jacobians re-run full forward solve; fresh coupler per solve |
| `optimize_layered` | **No** | Uses adjoint which assumes constant k/ρ/cp and no Q |

### Available Material YAML Models

Located in `MaterialStateSolver/tps_material_db/models/`:

| Material | Category | Key mechanism |
|----------|----------|---------------|
| `pica` | ablator | 2-stage phenolic pyrolysis, Bruggeman effective-medium k |
| `sla561v` | ablator | silicone pyrolysis + filler |
| `cork_p50` | ablator | cork + binder pyrolysis |
| `li900` | insulator | 93.5% porosity, radiation-surrogate k (T-dependent) |
| `li2200` | insulator | 84% porosity, radiation-surrogate k |
| `alumina` | ceramic | phonon scattering k(T), no reactions |
| `sic_sic_cmc` | ceramic | 10% porosity, effective medium |
| `carbon_carbon` | composite | Arrhenius carbon oxidation |
| `rcc` | composite | SiC coating → SiO₂, then carbon oxidation |
| `steel_304` | metal | T-dependent k and cp only, no reactions |
