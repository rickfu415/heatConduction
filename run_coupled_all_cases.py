"""
run_coupled_all_cases.py — run_all_cases with MaterialEngine coupling ON.

Same layered TPS cases as run_all_cases.py, but each hc.solve call receives
a MaterialEngine-backed material_hook so k/rho/cp/Q evolve with temperature
from kinetic YAML models in MaterialStateSolver/tps_material_db/models/.

Sections:
  A. Forward solve (1, 2, 3 layers) — coupled.
  B. optimize_layered (adjoint) — SKIPPED (adjoint assumes constant props).
  C. optimize_mass_slsqp (SLSQP + FD) — coupled.
  D. Over-designed → mass reduction — coupled.

Run:
    python run_coupled_all_cases.py
"""

import sys
from pathlib import Path

# sys.path hack — no install needed. Points at the sibling MaterialStateSolver repo.
_MSS_ROOT = Path(__file__).resolve().parent.parent.parent / 'MaterialStateSolver'
if _MSS_ROOT.exists() and str(_MSS_ROOT) not in sys.path:
    sys.path.insert(0, str(_MSS_ROOT))

import numpy as np
import parameter
import heatConduction as hc
import optimizer
from material_coupling import make_layered_coupler

# ── shared settings ─────────────────────────────────────────────────────────

T_BW_TARGET  = 450.0   # K
T_MIN        = 0.001   # m

STACKS = {
    1: {'materials': ['pica'],
        'thicknesses': np.array([0.050])},
    2: {'materials': ['pica', 'steel_304'],
        'thicknesses': np.array([0.035, 0.003])},
    3: {'materials': ['pica', 'li900', 'steel_304'],
        'thicknesses': np.array([0.030, 0.005, 0.001])},
}

SERVICE_TEMPS = {
    1: [np.inf],
    2: [np.inf, np.inf],
    3: [np.inf, 1530.0, 1200.0],
}


def make_para(n_layers):
    cfg = STACKS[n_layers]
    para = parameter.main()
    para['materials'] = list(cfg['materials'])
    t = cfg['thicknesses'].copy()
    para['layerThicknesses'] = t.copy()
    para['length'] = float(t.sum())
    para = parameter.load_materials(para)
    if n_layers == 1:
        props = parameter.load_material(cfg['materials'][0])
        para['material function'] = 'layered'
        para['numberOfLayers'] = 1
        para['layerConductivities'] = np.array([props['conductivity']])
        para['layerDensities']      = np.array([props['density']])
        para['layerHeatCapacities'] = np.array([props['heat_capacity']])
    para['back_wall_temperature_target'] = T_BW_TARGET
    parameter.normalize_conductivity(para)
    return para


def build_hook(para):
    """Fresh coupler → fresh material state per hc.solve call."""
    coupler = make_layered_coupler(para, list(para['materials']))
    return coupler.hook


def section(title):
    print('\n' + '='*70)
    print('  ' + title)
    print('='*70)


def fmt_mm(arr):
    return (np.asarray(arr) * 1000).round(2).tolist()


# ═══════════════════════════════════════════════════════════════════════════
# A. FORWARD SOLVE (coupled)
# ═══════════════════════════════════════════════════════════════════════════

for n in [1, 2, 3]:
    cfg = STACKS[n]
    mats = ' + '.join(cfg['materials'])
    section(f'A. FORWARD SOLVE (coupled) — {n} layer(s) ({mats})')
    para = make_para(n)
    TProfile, cache = hc.solve(para, verbose=False, material_hook=build_hook(para))
    T_bw = float(np.max(TProfile[-1, :]))
    print(f'  T_backwall_max = {T_bw:.1f} K  |  t = {fmt_mm(para["layerThicknesses"])} mm')


# ═══════════════════════════════════════════════════════════════════════════
# B. optimize_layered — SKIPPED (adjoint is invalid with evolving properties)
# ═══════════════════════════════════════════════════════════════════════════

section('B. optimize_layered — SKIPPED')
print('  Adjoint in differential.py assumes constant k/rho/cp. Use')
print('  optimize_mass_slsqp (FD-based) for coupled optimization.')


# ═══════════════════════════════════════════════════════════════════════════
# C. MINIMIZE MASS (SLSQP, coupled)
# ═══════════════════════════════════════════════════════════════════════════

for n in [1, 2, 3]:
    cfg = STACKS[n]
    mats = ' + '.join(cfg['materials'])
    svc = SERVICE_TEMPS[n]
    svc_str = ', '.join(
        f'T_{cfg["materials"][i]}≤{int(svc[i])} K' if not np.isinf(svc[i])
        else f'T_{cfg["materials"][i]}=free' for i in range(n)
    )
    section(f'C. MINIMIZE MASS (coupled) — {n} layer(s) ({mats})\n     '
            f'T_bw≤{int(T_BW_TARGET)} K  |  {svc_str}')

    para = make_para(n)
    para['material_engine_yamls'] = list(para['materials'])  # enables hook in optimizer
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
# D. OVER-DESIGNED → MASS REDUCTION (coupled)
# ═══════════════════════════════════════════════════════════════════════════

STACKS_OVER = {
    1: {'materials': ['pica'],
        'thicknesses': np.array([0.080])},
    2: {'materials': ['pica', 'steel_304'],
        'thicknesses': np.array([0.060, 0.005])},
    3: {'materials': ['pica', 'li900', 'steel_304'],
        'thicknesses': np.array([0.050, 0.010, 0.002])},
}

for n in [1, 2, 3]:
    cfg_over = STACKS_OVER[n]
    mats = ' + '.join(cfg_over['materials'])
    svc = SERVICE_TEMPS[n]
    section(f'D. OVER-DESIGNED → MASS REDUCTION (coupled) — {n} layer(s) ({mats})')

    para = make_para(n)
    t_over = cfg_over['thicknesses'].copy()
    para['layerThicknesses'] = t_over.copy()
    para['length'] = float(t_over.sum())
    parameter.normalize_conductivity(para)

    TProfile_init, _ = hc.solve(para, verbose=False, material_hook=build_hook(para))
    T_bw_init = float(np.max(TProfile_init[-1, :]))
    rho = np.atleast_1d(np.asarray(para['layerDensities'], dtype=float))
    mass_init = float(np.dot(rho, t_over))
    print(f'  Initial design:   t = {fmt_mm(t_over)} mm')
    print(f'                    mass = {mass_init:.3f} kg/m²  (static-density basis)')
    print(f'                    T_bw = {T_bw_init:.1f} K  (limit={int(T_BW_TARGET)} K)')

    para['material_engine_yamls'] = list(para['materials'])
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
          f'(saved {mass_saved_pct:.1f}%)')
    print(f'                    T_bw = {r["T_bw_max"]:.1f} K')
    print(f'                    T_layer_max = {r["T_layer_max"].round(1).tolist()} K')
    print(f'  Status:  {res.message}')
    print(f'  Constraints: {"ALL SATISFIED" if all_ok else "VIOLATED — check above"}')

print('\n' + '='*70)
print('  All coupled cases complete.')
print('='*70)
