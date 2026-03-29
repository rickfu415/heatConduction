# TPS Optimizer Guide

The optimizer exposes two entry points depending on the engineering question.

---

## At a Glance

| Question | Entry point | What it optimizes |
|---|---|---|
| "My TPS is too thin — what thickness meets the backwall T target?" | `optimize_layered` | Layer thicknesses via adjoint gradient descent |
| "My TPS is over-built — what's the minimum mass that still works?" | `optimize_mass_slsqp` | Layer thicknesses via SLSQP, three simultaneous constraints |

---

## `optimize_mass_slsqp` — Constrained Mass Minimization

### The Problem It Solves

```
minimize    m(t) = Σ ρᵢ · tᵢ            (areal mass, kg/m²)

subject to  T_backwall_max(t) ≤ T_bw_limit           [constraint 1]
            T_layer_i_max(t)  ≤ layer_service_temps[i] for each layer i  [constraint 2]
            tᵢ ≥ t_min[i]                             [constraint 3]
```

All three constraint types are active simultaneously. The solver finds the **thinnest stack
that keeps the back wall cold enough and no material above its melting point.**

### Signature

```python
res = optimizer.optimize_mass_slsqp(
    para_base,            # pd.Series — parameter set (never mutated)
    t0,                   # np.ndarray (n_layers,)  — initial thicknesses, metres
    T_bw_limit,           # float — back wall temperature limit, K
    layer_service_temps,  # list[float], length n_layers — per-layer limits, K
                          #   use np.inf to leave a layer unconstrained
    t_min,                # float or np.ndarray — per-layer lower bounds, m
    t_max=None,           # optional upper bounds (packaging constraint)
    tol=1e-6,             # SLSQP convergence tolerance
    max_iter=300,         # maximum iterations
    callback=None,        # called each iteration with progress dict
)
```

### Minimal Example — 2 Layers (pica + steel)

```python
import numpy as np
import parameter, optimizer

# --- build parameter set ---
para = parameter.main()
para['materials'] = ['pica', 'steel_304']
para['layerThicknesses'] = np.array([0.035, 0.003])
para['length'] = 0.038
para = parameter.load_materials(para)
para['back_wall_temperature_target'] = 450.0
parameter.normalize_conductivity(para)

# --- run optimizer ---
res = optimizer.optimize_mass_slsqp(
    para_base          = para,
    t0                 = np.array([0.035, 0.003]),   # initial design
    T_bw_limit         = 450.0,                      # K
    layer_service_temps= [np.inf, np.inf],            # no service-T limits here
    t_min              = np.array([0.001, 0.001]),    # 1 mm manufacturing floor
)

print(f"Optimized mass:  {res.fun:.3f} kg/m²")
print(f"Optimized layers: {(res.x * 1000).round(2)} mm")
```

### 3-Layer Example With Per-Layer Service Temperature Limits

```python
import numpy as np
import parameter, optimizer

# --- build parameter set ---
para = parameter.main()
para['materials'] = ['pica', 'li900', 'steel_304']
para['layerThicknesses'] = np.array([0.050, 0.010, 0.002])
para['length'] = 0.062
para = parameter.load_materials(para)
para['back_wall_temperature_target'] = 450.0
parameter.normalize_conductivity(para)

# --- constraints ---
T_bw_limit          = 450.0                         # K — structural back wall limit
layer_service_temps = [np.inf, 1530.0, 1200.0]      # LI-900 glass ≤ 1530 K, steel ≤ 1200 K
t_min               = np.array([0.001, 0.001, 0.001])

# --- run ---
res = optimizer.optimize_mass_slsqp(
    para_base          = para,
    t0                 = np.array([0.050, 0.010, 0.002]),
    T_bw_limit         = T_bw_limit,
    layer_service_temps= layer_service_temps,
    t_min              = t_min,
)

# --- verify ---
r = optimizer.compute_thermal_sensitivities(res.para_opt)
print(f"Optimized mass:      {res.fun:.3f} kg/m²")
print(f"t_opt (mm):          {(res.x * 1000).round(2).tolist()}")
print(f"T_backwall_max:      {r['T_bw_max']:.1f} K  (limit {T_bw_limit:.0f} K)")
print(f"T_layer_max (K):     {r['T_layer_max'].round(1).tolist()}")
print(f"Service T limits(K): {layer_service_temps}")
```

Expected output (approximate):
```
Optimized mass:      14.03 kg/m²
t_opt (mm):          [19.0, 7.6, 1.0]
T_backwall_max:      450.0 K  (limit 450 K)
T_layer_max (K):     [3142.4, 1530.0, 450.1]
Service T limits(K): [inf, 1530.0, 1200.0]
```

The optimizer simultaneously drove:
- **Back wall** to exactly the 450 K limit (constraint 1 active)
- **LI-900** to exactly its 1530 K service limit (constraint 2 active for layer 1)
- **Steel** well below 1200 K (constraint 2 inactive — back wall constraint dominates)
- **Each layer** above 1 mm (constraint 3 active for steel)

---

## Return Value

`res` is a `scipy.optimize.OptimizeResult` with extra attributes:

| Attribute | Type | Description |
|---|---|---|
| `res.x` | `np.ndarray (n_layers,)` | Optimized thicknesses in **metres** |
| `res.fun` | `float` | Optimized areal mass, kg/m² |
| `res.message` | `str` | SLSQP exit status |
| `res.success` | `bool` | True if converged |
| `res.para_opt` | `pd.Series` | Full parameter set at the optimum (ready for `hc.solve`) |
| `res.history` | `list[dict]` | Per-iteration log (see below) |

### Iteration History

Each entry in `res.history` is a dict:

```python
{
    'iter':                 int,
    'mass_kg_m2':           float,
    't_layers_mm':          list[float],
    'T_bw_max':             float,
    'T_layer_max':          list[float],
    'constraint_violations': {
        'bw':    float,           # max(0, T_bw_max - T_bw_limit)
        'layer': list[float],     # per constrained layer
    },
}
```

---

## `compute_thermal_sensitivities` — Inspect Any Design Point

Use this to evaluate temperatures and gradients at any parameter set, not just the optimum.

```python
r = optimizer.compute_thermal_sensitivities(para)

r['T_bw_max']       # float         — peak back wall T over the whole transient
r['T_layer_max']    # (n_layers,)   — peak T reached anywhere inside each layer
r['dT_bw_dt']       # (n_layers,)   — ∂T_backwall_max / ∂tᵢ  (K/m)
r['dT_layer_dt']    # (n_layers, n_layers) — [i,j] = ∂T_layer_i_max / ∂tⱼ  (K/m)
r['TProfile']       # (n_nodes, n_steps+1) — full temperature field
```

`dT_bw_dt[j] < 0` means increasing layer j thickness cools the back wall — the usual
case for an insulator. `dT_bw_dt[j] > 0` (rare) means the layer traps heat near the
back wall.

`dT_layer_dt[i, j]` captures cross-layer coupling: making the outer insulator (layer 0)
thicker shields all inner layers, so `dT_layer_dt[i, 0] < 0` for all i > 0.

---

## How the Gradients Are Computed

```
optimize_mass_slsqp
│
├── objective gradient   ∂m/∂tᵢ = ρᵢ          (analytic, trivial)
│
└── constraint Jacobian  ∂T/∂tⱼ               (finite differences)
      │
      ├── base forward solve:   T(t)
      └── for each layer j:
            perturb tⱼ → tⱼ + ε
            solve:  T(t + εeⱼ)
            FD:     ∂T/∂tⱼ ≈ [T(t + εeⱼ) − T(t)] / ε
```

**Why not the adjoint here?** The adjoint in `differential.py` computes
`dJ/d{k, ρ, cp}` per node at O(1) cost — it is used by `optimize_layered` via
`chain_rule`. The chain rule maps boundary-shift derivatives (zero-sum by
construction) to layer thicknesses, which is correct for redistribution but
wrong for independent thickness changes needed by SLSQP. Extending the adjoint
to independent thickness sensitivities requires differentiating through ghost
cells and non-uniform grid spacing at every layer interface — essentially
duplicating the full forward assembly. For 2–5 layers, n_layers extra forward
solves per SLSQP iteration are negligible.

---

## `optimize_layered` — Hit a Back Wall Temperature Target

Use when the design is **under-insulated** (T_backwall > target) and the goal is to
grow and redistribute layer thicknesses until the target is met.

```python
para_opt, history = optimizer.optimize_layered(
    para,         # parameter set — mutated in-place and also returned
    max_iter=200,
    tol=1.0,      # stop when loss = 0.5*(T_bw - target)² < tol
    callback=None,
)

# Read result
T_final = history[-1]['T_L']            # back wall T at last iteration
t_opt   = para_opt['layerThicknesses']  # optimized thicknesses
```

This optimizer uses the full adjoint pipeline:

```
forward solve → differential.main() → chain_rule → gradient descent
```

It does **not** minimize mass and does not enforce per-layer service temperature limits.
Use `optimize_mass_slsqp` for those requirements.

---

## Choosing the Right Optimizer

```
Is T_backwall > target?  (under-insulated, need to add material)
    YES → optimize_layered
    NO  → optimize_mass_slsqp

Do you have per-layer service temperature limits?
    YES → optimize_mass_slsqp  (only optimizer that enforces them)

Do you want minimum-mass output?
    YES → optimize_mass_slsqp
```
