"""
run_all_cases.py — comprehensive demonstration of all solver capabilities.

Cases:
  A. Forward solve only         (1, 2, 3 layers)
  B. optimize_layered           (1, 2, 3 layers)  — hit T_bw target, adjoint gradient descent
  C. optimize_mass_slsqp        (1, 2, 3 layers)  — minimize mass s.t. T_bw + service-T constraints

Run:
    python run_all_cases.py
"""

import numpy as np
import os
import parameter
import heatConduction as hc
import optimizer

# ── shared settings ─────────────────────────────────────────────────────────

T_BW_TARGET  = 450.0   # K  back-wall temperature target / limit
T_MIN        = 0.001   # m  manufacturing minimum thickness per layer (1 mm)

# Layer stacks for forward solve (A) and SLSQP (C)
STACKS = {
    1: {'materials': ['pica'],
        'thicknesses': np.array([0.050])},
    2: {'materials': ['pica', 'steel_304'],
        'thicknesses': np.array([0.035, 0.003])},
    3: {'materials': ['pica', 'li900', 'steel_304'],
        'thicknesses': np.array([0.030, 0.005, 0.001])},
}

# Starting stacks for optimize_layered (B): must start with T_bw > target so
# the growth/redistribution logic reduces T_bw toward the target.  Stacks A/C
# can start anywhere (auto-init handles SLSQP; forward solve is just diagnostic).
STACKS_B = {
    1: {'materials': ['pica'],
        'thicknesses': np.array([0.020])},          # thin → T_bw >> 450 K
    2: {'materials': ['pica', 'steel_304'],
        'thicknesses': np.array([0.010, 0.001])},   # thin pica → T_bw >> 450 K
    3: {'materials': ['pica', 'li900', 'steel_304'],
        'thicknesses': np.array([0.010, 0.002, 0.001])}, # under-insulated → T_bw >> 450 K
}

# Per-layer service temperature limits for SLSQP (inf = unconstrained)
SERVICE_TEMPS = {
    1: [np.inf],
    2: [np.inf, np.inf],
    3: [np.inf, 1530.0, 1200.0],
}


# ── helpers ──────────────────────────────────────────────────────────────────

def make_para(n_layers):
    """Build a parameter Series for an n-layer stack."""
    cfg = STACKS[n_layers]
    para = parameter.main()
    para['materials'] = list(cfg['materials'])
    t = cfg['thicknesses'].copy()
    para['layerThicknesses'] = t.copy()
    para['length'] = float(t.sum())
    para = parameter.load_materials(para)

    if n_layers == 1:
        # Upgrade 'constant' path to 'layered' so FD Jacobian works for SLSQP
        props = parameter.load_material(cfg['materials'][0])
        para['material function'] = 'layered'
        para['numberOfLayers'] = 1
        para['layerConductivities'] = np.array([props['conductivity']])
        para['layerDensities']      = np.array([props['density']])
        para['layerHeatCapacities'] = np.array([props['heat_capacity']])

    para['back_wall_temperature_target'] = T_BW_TARGET
    parameter.normalize_conductivity(para)
    return para


def section(title):
    print('\n' + '='*70)
    print('  ' + title)
    print('='*70)


def fmt_mm(arr):
    return (np.asarray(arr) * 1000).round(2).tolist()


# ═══════════════════════════════════════════════════════════════════════════
# A. FORWARD SOLVE ONLY
# ═══════════════════════════════════════════════════════════════════════════

section('A. FORWARD SOLVE — 1 layer (pica)')
para = make_para(1)
TProfile, cache = hc.solve(para)
T_bw = float(np.max(TProfile[-1, :]))
print(f'  T_backwall_max = {T_bw:.1f} K  |  t = {fmt_mm(para["layerThicknesses"])} mm')

section('A. FORWARD SOLVE — 2 layers (pica + steel_304)')
para = make_para(2)
TProfile, cache = hc.solve(para)
T_bw = float(np.max(TProfile[-1, :]))
print(f'  T_backwall_max = {T_bw:.1f} K  |  t = {fmt_mm(para["layerThicknesses"])} mm')

section('A. FORWARD SOLVE — 3 layers (pica + li900 + steel_304)')
para = make_para(3)
TProfile, cache = hc.solve(para)
T_bw = float(np.max(TProfile[-1, :]))
print(f'  T_backwall_max = {T_bw:.1f} K  |  t = {fmt_mm(para["layerThicknesses"])} mm')


# ═══════════════════════════════════════════════════════════════════════════
# B. OPTIMIZE BACKWALL T  (optimize_layered — adjoint gradient descent)
# ═══════════════════════════════════════════════════════════════════════════

for n in [1, 2, 3]:
    cfg = STACKS_B[n]
    mats = ' + '.join(cfg['materials'])
    section(f'B. OPTIMIZE T_bw → {T_BW_TARGET} K — {n} layer(s) ({mats})')
    para = make_para(n)
    # Override thicknesses with the under-insulated starting stack
    t_b = cfg['thicknesses'].copy()
    para['layerThicknesses'] = t_b.copy()
    para['length'] = float(t_b.sum())
    parameter.normalize_conductivity(para)
    para_opt, hist = optimizer.optimize_layered(para)
    T_final = hist[-1]['T_L']
    print(f'  T_bw_final = {T_final:.1f} K  |  t_opt = {fmt_mm(para_opt["layerThicknesses"])} mm')


# ═══════════════════════════════════════════════════════════════════════════
# C. MINIMIZE MASS  (optimize_mass_slsqp — SLSQP + FD Jacobian)
# ═══════════════════════════════════════════════════════════════════════════

for n in [1, 2, 3]:
    cfg = STACKS[n]
    mats = ' + '.join(cfg['materials'])
    svc = SERVICE_TEMPS[n]
    svc_str = ', '.join(
        f'T_{cfg["materials"][i]}≤{int(svc[i])} K' if not np.isinf(svc[i]) else f'T_{cfg["materials"][i]}=free'
        for i in range(n)
    )
    section(f'C. MINIMIZE MASS — {n} layer(s) ({mats})\n     T_bw≤{int(T_BW_TARGET)} K  |  {svc_str}')

    para = make_para(n)
    t0 = cfg['thicknesses'].copy()
    t_min_arr = np.full(n, T_MIN)

    res = optimizer.optimize_mass_slsqp(
        para_base=para,
        t0=t0,
        T_bw_limit=T_BW_TARGET,
        layer_service_temps=svc,
        t_min=t_min_arr,
    )

    r = optimizer.compute_thermal_sensitivities(res.para_opt)
    print(f'\n  Status:         {res.message}')
    print(f'  Optimized mass: {res.fun:.3f} kg/m²')
    print(f'  t_opt (mm):     {fmt_mm(res.x)}')
    print(f'  T_bw_max:       {r["T_bw_max"]:.1f} K  (limit={int(T_BW_TARGET)} K)')
    print(f'  T_layer_max:    {r["T_layer_max"].round(1).tolist()} K')
    viol_bw = max(0.0, r['T_bw_max'] - T_BW_TARGET)
    viol_svc = [max(0.0, r['T_layer_max'][i] - svc[i]) for i in range(n) if not np.isinf(svc[i])]
    all_ok = viol_bw < 2.0 and all(v < 2.0 for v in viol_svc)
    print(f'  Constraints:    {"ALL SATISFIED" if all_ok else "VIOLATED — check above"}')


# ═══════════════════════════════════════════════════════════════════════════
# D. OVER-DESIGNED → MASS REDUCTION  (explicit before/after demo)
#    Start from a thick, over-insulated design (T_bw << limit) and show
#    how much mass SLSQP recovers while still satisfying all constraints.
# ═══════════════════════════════════════════════════════════════════════════

# Over-designed configs — same as STACKS (high thickness, T_bw well below 450 K)
STACKS_OVER = {
    1: {'materials': ['pica'],
        'thicknesses': np.array([0.080])},          # 80 mm — very thick
    2: {'materials': ['pica', 'steel_304'],
        'thicknesses': np.array([0.060, 0.005])},   # 60 mm pica + 5 mm steel — heavy
    3: {'materials': ['pica', 'li900', 'steel_304'],
        'thicknesses': np.array([0.050, 0.010, 0.002])},  # 50+10+2 mm — over-built
}

for n in [1, 2, 3]:
    cfg_over = STACKS_OVER[n]
    mats = ' + '.join(cfg_over['materials'])
    svc = SERVICE_TEMPS[n]
    section(f'D. OVER-DESIGNED → MASS REDUCTION — {n} layer(s) ({mats})')

    para = make_para(n)
    t_over = cfg_over['thicknesses'].copy()
    para['layerThicknesses'] = t_over.copy()
    para['length'] = float(t_over.sum())
    parameter.normalize_conductivity(para)

    # ── Initial state ────────────────────────────────────────────────────
    TProfile_init, _ = hc.solve(para, verbose=False)
    T_bw_init = float(np.max(TProfile_init[-1, :]))
    rho = np.atleast_1d(np.asarray(para['layerDensities'], dtype=float))
    mass_init = float(np.dot(rho, t_over))
    print(f'  Initial design:   t = {fmt_mm(t_over)} mm')
    print(f'                    mass = {mass_init:.3f} kg/m²')
    print(f'                    T_bw = {T_bw_init:.1f} K  (limit={int(T_BW_TARGET)} K'
          f'  →  over-designed by {T_BW_TARGET - T_bw_init:.1f} K)')

    # ── SLSQP minimize mass ───────────────────────────────────────────────
    t_min_arr = np.full(n, T_MIN)
    res = optimizer.optimize_mass_slsqp(
        para_base=para,
        t0=t_over,
        T_bw_limit=T_BW_TARGET,
        layer_service_temps=svc,
        t_min=t_min_arr,
    )

    r = optimizer.compute_thermal_sensitivities(res.para_opt)
    mass_saved_pct = 100.0 * (mass_init - res.fun) / mass_init
    viol_bw = max(0.0, r['T_bw_max'] - T_BW_TARGET)
    viol_svc = [max(0.0, r['T_layer_max'][i] - svc[i]) for i in range(n) if not np.isinf(svc[i])]
    all_ok = viol_bw < 2.0 and all(v < 2.0 for v in viol_svc)

    print(f'\n  Optimized design: t = {fmt_mm(res.x)} mm')
    print(f'                    mass = {res.fun:.3f} kg/m²  '
          f'(saved {mass_saved_pct:.1f}% of initial mass)')
    print(f'                    T_bw = {r["T_bw_max"]:.1f} K')
    print(f'                    T_layer_max = {r["T_layer_max"].round(1).tolist()} K')
    print(f'  Status:  {res.message}')
    print(f'  Constraints: {"ALL SATISFIED" if all_ok else "VIOLATED — check above"}')

print('\n' + '='*70)
print('  All cases complete.')
print('='*70)
