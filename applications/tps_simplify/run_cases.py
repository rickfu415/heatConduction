"""
TPS single-layer comparison — evolving material properties via MaterialEngine.

Runs a 1D transient heat conduction solve for a single-layer slab of each
TPS material under identical boundary conditions:
  - hot face: constant applied heat flux + grey-body re-radiation
  - back face: adiabatic
  - duration, thickness, IC: all shared across materials

For each material we also run a static-JSON baseline (no state evolution)
for comparison. Results are serialized and plotted.

Run:
    cd .../heatConduction/applications/tps_simplify
    python run_cases.py

Requires MaterialStateSolver on PYTHONPATH (or pip-installed as material_engine).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Allow running this script from its folder by adding heatConduction/ root to sys.path.
_HC_ROOT = Path(__file__).resolve().parents[2]
if str(_HC_ROOT) not in sys.path:
    sys.path.insert(0, str(_HC_ROOT))

# Ensure MaterialStateSolver is importable.
_MSS_ROOT = _HC_ROOT.parent.parent / 'MaterialStateSolver'
if _MSS_ROOT.exists() and str(_MSS_ROOT) not in sys.path:
    sys.path.insert(0, str(_MSS_ROOT))

import parameter                             # noqa: E402
import heatConduction as hc                  # noqa: E402
from material_coupling import make_layered_coupler  # noqa: E402


# ── Test matrix ────────────────────────────────────────────────────────────

# Single-layer slabs under identical BCs. Thickness and flux are chosen so
# every material sees meaningful thermal response within the duration.
CASES = [
    {'name': 'pica',          'label': 'PICA (ablator)',          'category': 'ablator'},
    {'name': 'sla561v',       'label': 'SLA-561V (ablator)',      'category': 'ablator'},
    {'name': 'cork_p50',      'label': 'Cork P50 (ablator)',      'category': 'ablator'},
    {'name': 'li900',         'label': 'LI-900 (insulator)',      'category': 'insulator'},
    {'name': 'alumina',       'label': 'Alumina (ceramic)',       'category': 'ceramic'},
    {'name': 'steel_304',     'label': '304 Steel (metal)',       'category': 'metal'},
]

SHARED = {
    'thickness_m':     0.030,         # 30 mm slab
    'hot_flux_Wm2':    1.0e6,         # 1 MW/m^2 incident
    'duration_s':      120.0,
    'IC_K':            298.0,
    'nodesPerLayer':   41,
    'CFL':             3.0,
    'emissivity':      0.9,
    'ambient_K':       298.0,
}

# Milder, longer study for state-history plots. At 1 MW/m² ablator pyrolysis
# completes in the first timestep (T_hot hits ~2100 K immediately), which
# produces uninformative history curves. A lower flux lets the reaction
# front propagate gradually so we can see component densities, porosity
# and k/ρ evolve in time.
HISTORY_SHARED = {
    'thickness_m':     0.030,
    'hot_flux_Wm2':    2.0e5,         # 200 kW/m^2 — mild re-entry soak
    'duration_s':      600.0,
    'IC_K':            298.0,
    'nodesPerLayer':   41,
    'CFL':             3.0,
    'emissivity':      0.9,
    'ambient_K':       298.0,
}

OUT_DIR = Path(__file__).parent / 'figures'
OUT_DIR.mkdir(exist_ok=True)


# ── Helpers ────────────────────────────────────────────────────────────────

def build_para(material_name, *, thickness_m, duration_s, nodesPerLayer,
               hot_flux_Wm2, IC_K, CFL, emissivity, ambient_K):
    """Build a 1-layer para Series from a material name (static JSON path)."""
    para = parameter.main()
    para['materials'] = [material_name]
    t = np.array([thickness_m])
    para['layerThicknesses'] = t.copy()
    para['length'] = float(t.sum())
    para = parameter.load_materials(para)

    # Upgrade single-layer 'constant' flow to 'layered' so material_coupling
    # can index per-layer arrays uniformly.
    props = parameter.load_material(material_name)
    para['material function'] = 'layered'
    para['numberOfLayers'] = 1
    para['layerConductivities'] = np.array([props['conductivity']])
    para['layerDensities']      = np.array([props['density']])
    para['layerHeatCapacities'] = np.array([props['heat_capacity']])
    para['nodesPerLayer'] = nodesPerLayer

    # BCs: hot-face flux + re-radiate, back-face adiabatic.
    para['x=0 type'] = 'heatFlux'
    para['x=0 value'] = float(hot_flux_Wm2)
    para['x=L type'] = 'heatFlux'
    para['x=L value'] = 0.0
    para['re-radiate'] = True
    para['emissivity'] = emissivity
    para['ambientTemperature'] = ambient_K
    para['IC value'] = float(IC_K)
    para['duration'] = float(duration_s)
    para['CFL'] = float(CFL)

    parameter.normalize_conductivity(para)
    return para


def run_one(material_name, *, coupled, shared):
    """Run forward solve for a single material, return a result dict."""
    para = build_para(material_name, **{k: shared[k] for k in (
        'thickness_m', 'duration_s', 'nodesPerLayer', 'hot_flux_Wm2',
        'IC_K', 'CFL', 'emissivity', 'ambient_K')})

    coupler = None
    hook = None
    if coupled:
        coupler = make_layered_coupler(para, [material_name])
        hook = coupler.hook

    TProfile, cache = hc.solve(para, verbose=False, material_hook=hook)

    dt = float(para['deltaTime'])
    nsteps = int(para['numberOfTimeStep'])
    t_arr = np.arange(nsteps + 1) * dt
    x_arr = np.concatenate(([0.0], np.cumsum(para['dx_array'])))

    return {
        'material': material_name,
        'coupled': coupled,
        't': t_arr,
        'x': x_arr,
        'T': TProfile,                     # shape (n_nodes, n_steps+1)
        'T_hot':  TProfile[0, :],
        'T_back': TProfile[-1, :],
        'para':   para,
        'coupler': coupler,                # None when static
    }


def snapshot_state(coupler, t_indices, TProfile):
    """For a coupled run, build a time-indexed snapshot of component densities
    and structural variables at selected hot/back nodes. Only makes sense for
    ablators; for inert materials most components don't change."""
    # The coupler holds only the FINAL state (state.n_points matches npl).
    # For intermediate snapshots we'd need to re-run stepwise; for this demo
    # we extract only the final state. For time-series of state we use a
    # separate stepping path in run_with_history().
    return {
        'components': {k: v.copy() for k, v in coupler.states[0].component_densities.items()},
        'structural': {k: v.copy() for k, v in coupler.states[0].structural_variables.items()},
    }


def run_with_history(material_name, shared):
    """Coupled run that records per-timestep state at hot and back nodes for
    plotting the pyrolysis history. Reimplements the hc.solve time loop.
    """
    para = build_para(material_name, **{k: shared[k] for k in (
        'thickness_m', 'duration_s', 'nodesPerLayer', 'hot_flux_Wm2',
        'IC_K', 'CFL', 'emissivity', 'ambient_K')})
    coupler = make_layered_coupler(para, [material_name])

    # Walk hc.solve's internals manually so we can snapshot state.
    para = parameter.normalize_conductivity(para)
    cache = hc.initialize(para)
    nsteps = int(para['numberOfTimeStep'])
    dt = float(para['deltaTime'])

    # Track 3 probe nodes: hot face, 25% depth, 50% depth.
    # (Back face rarely reacts — we know that; instead capture the
    # propagating pyrolysis front at shallow/mid depth.)
    npl = int(para['nodesPerLayer'])
    PROBES = [0, npl // 4, npl // 2]
    PROBE_LABELS = ['hot face (x=0)', '25% depth', '50% depth']

    comp_names = list(coupler.states[0].component_densities.keys())
    history = {c: np.zeros((nsteps + 1, len(PROBES))) for c in comp_names}
    struct_names = list(coupler.states[0].structural_variables.keys())
    structs = {c: np.zeros((nsteps + 1, len(PROBES))) for c in struct_names}
    k_hist = np.zeros((nsteps + 1, len(PROBES)))
    rho_hist = np.zeros((nsteps + 1, len(PROBES)))
    Q_hist = np.zeros((nsteps + 1, len(PROBES)))

    def snap(ts):
        for c in comp_names:
            arr = coupler.states[0].component_densities[c]
            for j, p in enumerate(PROBES):
                history[c][ts, j] = arr[p]
        for c in struct_names:
            arr = coupler.states[0].structural_variables[c]
            for j, p in enumerate(PROBES):
                structs[c][ts, j] = arr[p]
        for j, p in enumerate(PROBES):
            k_hist[ts, j]   = para['conductivity'][p]
            rho_hist[ts, j] = para['density'][p]
        Q = para.get('volumetricHeatSource', None)
        if Q is not None:
            for j, p in enumerate(PROBES):
                Q_hist[ts, j] = Q[p]

    snap(0)
    for ts in range(1, nsteps + 1):
        cache['ts'] = ts
        cache = hc.newtonIteration(para, cache, verbose=False)
        coupler.hook(para, cache, ts)
        cache = hc.storeUpdateResult(cache)
        snap(ts)

    t_arr = np.arange(nsteps + 1) * dt
    return {
        'material': material_name,
        't': t_arr,
        'components': history,
        'structural': structs,
        'k_hist': k_hist, 'rho_hist': rho_hist, 'Q_hist': Q_hist,
        'TProfile': cache['TProfile'],
        'probe_labels': PROBE_LABELS,
    }


# ── Driver ─────────────────────────────────────────────────────────────────

def main():
    print('='*70)
    print(' TPS single-layer comparison — evolving vs static materials')
    print(' slab {:.1f} mm | flux {:.2f} MW/m² + re-radiate | {:.0f} s'.format(
          SHARED['thickness_m']*1000, SHARED['hot_flux_Wm2']/1e6, SHARED['duration_s']))
    print('='*70)

    results_coupled = []
    results_static = []
    for case in CASES:
        name = case['name']
        print(f'  running {case["label"]:<28} ', end='', flush=True)
        r_cpl = run_one(name, coupled=True,  shared=SHARED)
        r_stc = run_one(name, coupled=False, shared=SHARED)
        r_cpl['label'] = case['label']; r_cpl['category'] = case['category']
        r_stc['label'] = case['label']; r_stc['category'] = case['category']
        results_coupled.append(r_cpl)
        results_static.append(r_stc)
        print(f'T_hot_max = {r_cpl["T_hot"].max():6.0f} K | '
              f'T_back_max = {r_cpl["T_back"].max():5.0f} K  '
              f'(static baseline back: {r_stc["T_back"].max():5.0f} K)')

    # Pyrolysis-state history — only meaningful for ablators
    histories = {}
    for case in CASES:
        if case['category'] == 'ablator':
            print(f'  recording state history for {case["label"]}... ', end='', flush=True)
            histories[case['name']] = run_with_history(case['name'], HISTORY_SHARED)
            print('done')

    # ── Plots ────────────────────────────────────────────────────────────
    plot_backwall(results_coupled, results_static)
    plot_hotface(results_coupled, results_static)
    plot_temperature_profiles(results_coupled)
    plot_pyrolysis_history(histories)
    plot_property_evolution(histories)

    # ── Summary table ────────────────────────────────────────────────────
    print()
    print('='*80)
    summary = pd.DataFrame([{
        'material': r['label'],
        'T_hot_max_coupled_K':  float(r['T_hot'].max()),
        'T_back_max_coupled_K': float(r['T_back'].max()),
        'T_back_max_static_K':  float(s['T_back'].max()),
        'delta_back_K':         float(r['T_back'].max() - s['T_back'].max()),
    } for r, s in zip(results_coupled, results_static)])
    print(summary.to_string(index=False, float_format=lambda x: f'{x:8.1f}'))
    summary.to_csv(OUT_DIR / 'summary.csv', index=False)
    print()
    print(f' Plots and summary.csv saved in: {OUT_DIR}')


# ── Plot routines ──────────────────────────────────────────────────────────

def _material_color(cat):
    return {
        'ablator':   '#c0392b',
        'insulator': '#2980b9',
        'ceramic':   '#8e44ad',
        'composite': '#16a085',
        'metal':     '#7f8c8d',
    }.get(cat, '#333333')


def plot_backwall(results_coupled, results_static):
    fig, ax = plt.subplots(figsize=(10, 6))
    for r_cpl, r_stc in zip(results_coupled, results_static):
        col = _material_color(r_cpl['category'])
        ax.plot(r_cpl['t'], r_cpl['T_back'], color=col, lw=2.2,
                label=f"{r_cpl['label']}  (coupled)")
        ax.plot(r_stc['t'], r_stc['T_back'], color=col, lw=1.2,
                ls='--', alpha=0.65, label=f"{r_stc['label']}  (static)")
    ax.set_xlabel('time (s)')
    ax.set_ylabel('back-face temperature (K)')
    ax.set_title('Back-face temperature: evolving vs static material properties\n'
                 f"slab {SHARED['thickness_m']*1000:.0f} mm | "
                 f"q_in={SHARED['hot_flux_Wm2']/1e6:.1f} MW/m² + re-radiate (ε={SHARED['emissivity']})")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=1, bbox_to_anchor=(1.02, 1.0), loc='upper left')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig_backwall_vs_time.png', dpi=150)
    plt.close(fig)


def plot_hotface(results_coupled, results_static):
    fig, ax = plt.subplots(figsize=(10, 6))
    for r_cpl, r_stc in zip(results_coupled, results_static):
        col = _material_color(r_cpl['category'])
        ax.plot(r_cpl['t'], r_cpl['T_hot'], color=col, lw=2.2,
                label=f"{r_cpl['label']}  (coupled)")
        ax.plot(r_stc['t'], r_stc['T_hot'], color=col, lw=1.2,
                ls='--', alpha=0.65)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('hot-face temperature (K)')
    ax.set_title('Hot-face temperature: evolving (solid) vs static (dashed)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=1, bbox_to_anchor=(1.02, 1.0), loc='upper left')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig_hotface_vs_time.png', dpi=150)
    plt.close(fig)


def plot_temperature_profiles(results_coupled):
    n = len(results_coupled)
    ncol = 3; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 3.2 * nrow),
                              sharex=True)
    axes = np.atleast_1d(axes).ravel()
    # Snapshots at 5 uniformly-spaced times
    for ax, r in zip(axes, results_coupled):
        col = _material_color(r['category'])
        nsteps = r['T'].shape[1] - 1
        idxs = np.linspace(0, nsteps, 6).astype(int)
        cmap = plt.cm.plasma(np.linspace(0.15, 0.90, len(idxs)))
        for i, ti in enumerate(idxs):
            ax.plot(r['x'] * 1000, r['T'][:, ti], color=cmap[i], lw=1.6,
                    label=f"t={r['t'][ti]:.0f} s")
        ax.set_title(r['label'], color=col, fontsize=10)
        ax.set_xlabel('x (mm)'); ax.set_ylabel('T (K)')
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc='upper right')
    # hide extras
    for ax in axes[n:]:
        ax.axis('off')
    fig.suptitle('In-depth temperature profiles (coupled material engine)')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig_temperature_profiles.png', dpi=150)
    plt.close(fig)


def plot_pyrolysis_history(histories):
    """Component density evolution at 3 depths for ablators (shows the
    propagating pyrolysis front)."""
    if not histories:
        return
    n = len(histories)
    any_h = next(iter(histories.values()))
    nprobes = len(any_h['probe_labels'])
    fig, axes = plt.subplots(n, nprobes, figsize=(4.0 * nprobes, 3.0 * n), sharex=True)
    axes = np.atleast_2d(axes)
    for row, (name, h) in enumerate(histories.items()):
        for col_idx, loc in enumerate(h['probe_labels']):
            ax = axes[row, col_idx]
            for comp, hist in h['components'].items():
                col = hist[:, col_idx]
                if col.max() < 1e-6 and col.min() < 1e-6:
                    continue
                ax.plot(h['t'], col, lw=1.6, label=comp)
            ax.set_title(f'{name} — {loc}', fontsize=10)
            if col_idx == 0:
                ax.set_ylabel('ρ (kg/m³)')
            ax.grid(alpha=0.3)
            if row == 0 and col_idx == nprobes - 1:
                ax.legend(fontsize=7, loc='upper right')
            if row == n - 1:
                ax.set_xlabel('t (s)')
    fig.suptitle('Ablator component density at 3 depths — milder soak '
                 f'({HISTORY_SHARED["hot_flux_Wm2"]/1e3:.0f} kW/m², '
                 f'{HISTORY_SHARED["duration_s"]:.0f} s)')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig_pyrolysis_history.png', dpi=150)
    plt.close(fig)


def plot_property_evolution(histories):
    """Bulk k and rho at 3 depths vs time for ablators."""
    if not histories:
        return
    n = len(histories)
    probe_colors = ['#c0392b', '#e67e22', '#2c3e50']
    fig, axes = plt.subplots(n, 2, figsize=(10, 3.0 * n), sharex=True)
    axes = np.atleast_2d(axes)
    for row, (name, h) in enumerate(histories.items()):
        for j, lbl in enumerate(h['probe_labels']):
            axes[row, 0].plot(h['t'], h['k_hist'][:, j], lw=1.8, label=lbl,
                              color=probe_colors[j % len(probe_colors)])
            axes[row, 1].plot(h['t'], h['rho_hist'][:, j], lw=1.8, label=lbl,
                              color=probe_colors[j % len(probe_colors)])
        axes[row, 0].set_title(f'{name} — conductivity k', fontsize=10)
        axes[row, 0].set_ylabel('k (W/m·K)')
        axes[row, 0].grid(alpha=0.3)
        axes[row, 1].set_title(f'{name} — bulk density ρ', fontsize=10)
        axes[row, 1].set_ylabel('ρ (kg/m³)')
        axes[row, 1].grid(alpha=0.3)
        if row == 0:
            axes[row, 0].legend(fontsize=7, loc='best')
            axes[row, 1].legend(fontsize=7, loc='best')
        if row == n - 1:
            axes[row, 0].set_xlabel('t (s)'); axes[row, 1].set_xlabel('t (s)')
    fig.suptitle('Bulk property evolution at 3 depths — milder soak '
                 f'({HISTORY_SHARED["hot_flux_Wm2"]/1e3:.0f} kW/m², '
                 f'{HISTORY_SHARED["duration_s"]:.0f} s)')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'fig_property_evolution.png', dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    main()
