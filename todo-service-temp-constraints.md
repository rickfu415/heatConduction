# TODO: Solver Optimization Capability Roadmap

## The Core Problem with the Current Approach

The current adjoint optimizer solves:

```
minimize    (T_backwall - T_target)²
```

This is the **wrong objective**. Mass grows as a side effect of hitting the target, with no
incentive to be efficient. The symptoms:

- "Design is heavy" → the optimizer doesn't care about mass, only about hitting the T target
- "Layer X is over-designed" → the optimizer redistributes thickness without penalizing excess
- Service temperature limits are ignored entirely

The engineering goal is actually a **constrained NLP**:

```
minimize    m(t) = Σ ρᵢ · tᵢ                    ← areal mass (kg/m²)

subject to  T_backwall_max(t) ≤ T_bw_limit       ← thermal performance
            T_layer_i_max(t)  ≤ T_service_i      ← material survivability (per layer)
            tᵢ ≥ tᵢ_min                          ← manufacturing lower bound
            t_total ≤ t_total_max (optional)      ← packaging constraint
```

where `T_backwall_max(t)` and `T_layer_i_max(t)` are evaluated by running the 1D transient
heat conduction solver.

---

## What the Solver Should Be Capable Of

### Capability 1: Mass-Minimization with Thermal Constraints (primary goal)

Reframe the core optimization as the constrained NLP above.

**Objective**: minimize total areal mass.
**Hard constraints**:
- `T_backwall_max ≤ T_bw_target` — back wall must stay below structural limit
- `T_layer_i_max ≤ T_service_i` for each layer — materials must not be destroyed

**Implementation approach — SLSQP via `scipy.optimize.minimize`**:
```python
from scipy.optimize import minimize, NonlinearConstraint, Bounds

def optimize_mass_slsqp(
    para_base: pd.Series,
    t0: np.ndarray,                    # initial thicknesses (m)
    T_bw_limit: float,                 # back wall temperature limit (K)
    layer_service_temps: list[float],  # per-layer service limits (K)
    t_min: np.ndarray,                 # lower bounds on each thickness (m)
    t_max: np.ndarray | None = None,   # optional upper bounds
    tol: float = 0.5,                  # constraint satisfaction tolerance (K)
    callback=None,
) -> OptimizeResult:
```

This is already in scipy — no new solver code needed. The key is providing:
1. **Objective function**: `f(t) = Σ ρᵢ tᵢ`
2. **Gradient of objective**: `df/dtᵢ = ρᵢ` (analytic — trivial)
3. **Constraint functions**: `g_bw(t) = T_bw_limit - T_backwall_max(t) ≥ 0`
   and `g_layer_i(t) = T_service_i - T_layer_i_max(t) ≥ 0`
4. **Constraint Jacobian**: `dg/dtⱼ = -dT/dtⱼ` (via adjoint — already partially implemented)

SLSQP handles all inequality constraints natively and converges in O(10–100) forward solves
for 2–5 layers.

---

### Capability 2: Gradient Computation for All Thermal Quantities

Currently the adjoint computes `dT_backwall/dtⱼ` only. Needed additions:

```python
def compute_thermal_sensitivities(
    para: pd.Series,
) -> dict:
    """
    Returns:
        dT_bw_dt   : np.ndarray shape (n_layers,)  — sensitivity of peak back wall T
        dT_layer_dt: np.ndarray shape (n_layers, n_layers) — dT_layer_i_max / dt_j
    """
```

`dT_layer_i_max / dt_j` is the sensitivity of layer i's peak temperature to the thickness
of layer j. This can be computed via:
- **Adjoint (Option A, preferred)**: extend the existing adjoint to track interior
  temperature sensitivities, not just the back wall. Accurate, O(1) extra cost.
- **Finite differences (Option B, fallback)**: perturb each `tⱼ`, re-run solver,
  compute difference. O(n_layers) extra solves per iteration — fine for 2–3 layers.

---

### Capability 3: Gradient-Free Search (for material selection or non-smooth problems)

When the decision space includes discrete material selection, gradient methods break down.
Add a gradient-free optimizer as an alternative engine:

```python
def optimize_mass_de(
    para_base: pd.Series,
    material_library: list[MaterialSpec],  # candidate materials per layer slot
    T_bw_limit: float,
    layer_service_temps: list[float],
    ...
) -> OptimizeResult:
```

**Suggested method**: Differential Evolution (`scipy.optimize.differential_evolution`)
- Handles mixed continuous (thickness) + discrete (material) spaces
- No gradients needed
- Naturally handles non-smooth, multi-modal landscapes
- Cost: O(population × generations) forward solves — typically 200–2000 evaluations

For thickness-only problems with fixed materials, DE is overkill. Reserve for material
selection or when the landscape is suspected non-smooth (e.g., ablating materials,
phase-change effects in future).

---

### Capability 4: Pareto Front for Multi-Objective Problems

When there are two competing objectives (e.g., minimize mass AND minimize back wall T),
expose a Pareto front sweep:

```python
def pareto_sweep(
    para_base: pd.Series,
    T_bw_targets: list[float],   # sweep over different back wall T limits
    ...
) -> list[OptimizeResult]:
```

Mechanically: run Capability 1 (SLSQP) at each `T_bw_target` value, collect (mass, T_bw)
pairs. This gives the Pareto front with minimal new code.

---

## Recommended Implementation Order

| Step | What | Method | New code in solver | Status |
|------|------|--------|--------------------|--------|
| 1 | Mass minimization, back wall T constraint only | SLSQP + FD Jacobian | `optimize_mass_slsqp()` entry point | ✅ Done (2026-03-28) |
| 2 | Add per-layer service T constraints | FD Jacobian (same call) | `compute_thermal_sensitivities()` | ✅ Done (2026-03-28) |
| 3 | Extend adjoint for interior T sensitivities | Adjoint extension | Would replace FD in Steps 1–2 | Deferred — see note |
| 4 | Material selection | Differential Evolution | `optimize_mass_de()` | Not started |
| 5 | Pareto front | SLSQP sweep | `pareto_sweep()` utility | Not started |

Steps 1–2 are complete. Steps 4–5 are the remaining high-value items.

**Note on Step 3:** The `chain_rule` function computes zero-sum boundary-shift derivatives
(correct for `optimize_layered` redistribution). Extending it analytically to independent
thickness variations for SLSQP requires differentiating through ghost cells, BC source terms,
and non-uniform dx at layer interfaces — essentially duplicating the entire forward assembly.
FD on the forward solver gives exact independent sensitivities at the cost of n_layers extra
forward solves per SLSQP iteration, which is negligible for 2–5 layers. Step 3 is deferred
unless solve speed becomes a bottleneck at larger layer counts.

---

## Interface Contract with the Web App

Once the solver is updated, the web app backend (`backend/solver/optimize.py`) needs to:

1. Expose a new endpoint concept: **mass minimization** (vs. current "hit T target")
2. The optimizer request should carry:
   - `objective: "minimize_mass"` (new) or `"hit_target"` (current)
   - `T_bw_limit: float` — back wall temperature limit
   - `layer_service_temps: list[float]` — per-layer service limits
   - `t_min_mm: list[float]` — per-layer minimum thickness (e.g., structural floor)
3. The SSE progress events should include:
   - `mass_kg_m2` — current areal mass at this iteration
   - `constraint_violations` — which constraints are violated and by how much
   - `layer_max_temps` — per-layer peak T at current thicknesses
4. The "done" event should include:
   - `layer_max_temps` — final per-layer peak T (from last forward solve)
   - `constraint_satisfied` — boolean per constraint
   - `mass_kg_m2_original` and `mass_kg_m2_optimized` — for the UI delta display

---

## Open Questions

1. **Minimum thickness bounds**: What is a physically meaningful `t_min` per layer?
   (e.g., 1 mm structural floor, or material-dependent)
2. **Starting point sensitivity**: SLSQP is gradient-based and can get stuck. Does the
   current design act as a good initial guess? Should we add a multi-start strategy?
3. **Constraint tolerance**: How tight should the back wall T constraint be enforced?
   `T_bw ≤ limit + 2K` is probably fine, but SLSQP needs a feasible starting point or
   good constraint scaling.
4. **Adjoint vs FD for layer sensitivities**: The existing adjoint computes back-wall
   sensitivity. Extending it to interior nodes is the right long-term choice but requires
   understanding the adjoint derivation deeply. FD is a safe first step.
