"""
Performance benchmark: heatConduction with vs without MaterialEngine coupling.

Measures wall-clock time for forward solves and SLSQP optimization across
materials, layer counts, and grid resolutions. Outputs a summary table and
a bar-chart comparison figure.

Run:
    cd .../heatConduction/applications/tps_simplify
    python run_benchmark.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

_HC_ROOT = Path(__file__).resolve().parents[2]
if str(_HC_ROOT) not in sys.path:
    sys.path.insert(0, str(_HC_ROOT))
_MSS_ROOT = _HC_ROOT.parent.parent / 'MaterialStateSolver'
if _MSS_ROOT.exists() and str(_MSS_ROOT) not in sys.path:
    sys.path.insert(0, str(_MSS_ROOT))

import parameter
import heatConduction as hc
import optimizer
from material_coupling import make_layered_coupler

OUT_DIR = Path(__file__).parent / 'figures'
OUT_DIR.mkdir(exist_ok=True)


# ── Helpers ────────────────────────────────────────────────────────────────

def build_para(materials, thicknesses, *, nodesPerLayer=21, duration=120.0,
               flux=1e6, CFL=3.0):
    """Build a layered para Series."""
    para = parameter.main()
    para['materials'] = list(materials)
    t = np.array(thicknesses, dtype=float)
    para['layerThicknesses'] = t.copy()
    para['length'] = float(t.sum())
    para = parameter.load_materials(para)

    n = len(materials)
    if n == 1:
        props = parameter.load_material(materials[0])
        para['material function'] = 'layered'
        para['numberOfLayers'] = 1
        para['layerConductivities'] = np.array([props['conductivity']])
        para['layerDensities']      = np.array([props['density']])
        para['layerHeatCapacities'] = np.array([props['heat_capacity']])

    para['nodesPerLayer'] = nodesPerLayer
    para['x=0 type'] = 'heatFlux'
    para['x=0 value'] = float(flux)
    para['x=L type'] = 'heatFlux'
    para['x=L value'] = 0.0
    para['re-radiate'] = True
    para['emissivity'] = 0.9
    para['ambientTemperature'] = 298.0
    para['IC value'] = 298.0
    para['duration'] = float(duration)
    para['CFL'] = float(CFL)
    para['back_wall_temperature_target'] = 450.0
    parameter.normalize_conductivity(para)
    return para


def time_forward(para, coupled):
    """Time a single forward solve. Returns (wall_seconds, n_timesteps, n_nodes)."""
    hook = None
    if coupled:
        coupler = make_layered_coupler(para, list(para['materials']))
        hook = coupler.hook

    t0 = time.perf_counter()
    TProfile, cache = hc.solve(para, verbose=False, material_hook=hook)
    elapsed = time.perf_counter() - t0

    nsteps = int(para['numberOfTimeStep'])
    nnodes = int(para['numberOfNode'])
    return elapsed, nsteps, nnodes


def time_slsqp(para, coupled, max_iter=5):
    """Time SLSQP with a fixed small number of iterations. Returns wall seconds."""
    n = len(np.atleast_1d(para['layerConductivities']))
    t0_arr = np.atleast_1d(np.asarray(para['layerThicknesses'], dtype=float))

    para_run = para.copy()
    if coupled:
        para_run['material_engine_yamls'] = list(para['materials'])

    t0 = time.perf_counter()
    optimizer.optimize_mass_slsqp(
        para_base=para_run,
        t0=t0_arr,
        T_bw_limit=450.0,
        layer_service_temps=[np.inf] * n,
        t_min=np.full(n, 0.001),
        max_iter=max_iter,
    )
    elapsed = time.perf_counter() - t0
    return elapsed


# ── Benchmark configurations ──────────────────────────────────────────────

FORWARD_CASES = [
    # (label, materials, thicknesses, nodesPerLayer)
    ('PICA 1L npl=21',        ['pica'],                      [0.050],               21),
    ('PICA 1L npl=41',        ['pica'],                      [0.050],               41),
    ('PICA 1L npl=81',        ['pica'],                      [0.050],               81),
    ('SLA-561V 1L npl=21',    ['sla561v'],                   [0.030],               21),
    ('LI-900 1L npl=21',      ['li900'],                     [0.030],               21),
    ('Steel 1L npl=21',       ['steel_304'],                 [0.005],               21),
    ('PICA+Steel 2L npl=21',  ['pica', 'steel_304'],         [0.035, 0.003],        21),
    ('3L stack npl=21',       ['pica', 'li900', 'steel_304'],[0.030, 0.005, 0.001], 21),
]

SLSQP_CASES = [
    # (label, materials, thicknesses, nodesPerLayer, max_iter)
    ('PICA 1L',       ['pica'],                      [0.050],               21, 3),
    ('PICA+Steel 2L', ['pica', 'steel_304'],         [0.035, 0.003],        21, 3),
    ('3L stack',      ['pica', 'li900', 'steel_304'],[0.030, 0.005, 0.001], 21, 3),
]


# ── Run benchmarks ────────────────────────────────────────────────────────

def main():
    print('='*80)
    print(' Performance benchmark: static vs MaterialEngine-coupled')
    print('='*80)

    # ── Forward solves ────────────────────────────────────────────────────
    print('\n--- Forward Solve ---')
    print(f'{"Case":<26} {"nodes":>5} {"steps":>5}  '
          f'{"static (s)":>10} {"coupled (s)":>11} {"overhead":>9}')
    print('-'*80)

    fwd_results = []
    for label, mats, thick, npl in FORWARD_CASES:
        para = build_para(mats, thick, nodesPerLayer=npl)

        t_static, nsteps, nnodes = time_forward(para.copy(), coupled=False)
        t_coupled, _, _          = time_forward(para.copy(), coupled=True)
        overhead = t_coupled / t_static if t_static > 0 else float('inf')

        fwd_results.append({
            'label': label, 'nodes': nnodes, 'steps': nsteps,
            't_static': t_static, 't_coupled': t_coupled, 'overhead': overhead,
        })
        print(f'{label:<26} {nnodes:>5} {nsteps:>5}  '
              f'{t_static:>10.3f} {t_coupled:>11.3f} {overhead:>8.1f}x')

    # ── SLSQP ────────────────────────────────────────────────────────────
    print('\n--- SLSQP Optimizer (fixed iterations) ---')
    print(f'{"Case":<26} {"iters":>5}  '
          f'{"static (s)":>10} {"coupled (s)":>11} {"overhead":>9}')
    print('-'*80)

    slsqp_results = []
    for label, mats, thick, npl, maxiter in SLSQP_CASES:
        para = build_para(mats, thick, nodesPerLayer=npl)

        t_static  = time_slsqp(para.copy(), coupled=False,  max_iter=maxiter)
        t_coupled = time_slsqp(para.copy(), coupled=True,   max_iter=maxiter)
        overhead = t_coupled / t_static if t_static > 0 else float('inf')

        slsqp_results.append({
            'label': label, 'iters': maxiter,
            't_static': t_static, 't_coupled': t_coupled, 'overhead': overhead,
        })
        print(f'{label:<26} {maxiter:>5}  '
              f'{t_static:>10.3f} {t_coupled:>11.3f} {overhead:>8.1f}x')

    # ── Summary ──────────────────────────────────────────────────────────
    overheads = [r['overhead'] for r in fwd_results]
    print(f'\nForward-solve overhead range: {min(overheads):.1f}x – {max(overheads):.1f}x')
    overheads_slsqp = [r['overhead'] for r in slsqp_results]
    print(f'SLSQP overhead range:        {min(overheads_slsqp):.1f}x – {max(overheads_slsqp):.1f}x')

    # ── Plots ────────────────────────────────────────────────────────────
    plot_forward_bar(fwd_results)
    plot_slsqp_bar(slsqp_results)
    plot_scaling(fwd_results)

    print(f'\nFigures saved in: {OUT_DIR}')


# ── Plot routines ─────────────────────────────────────────────────────────

def plot_forward_bar(results):
    labels = [r['label'] for r in results]
    t_s = [r['t_static'] for r in results]
    t_c = [r['t_coupled'] for r in results]

    x = np.arange(len(labels))
    w = 0.35

    fig, ax = plt.subplots(figsize=(11, 5))
    bars_s = ax.bar(x - w/2, t_s, w, label='static (JSON)', color='#95a5a6')
    bars_c = ax.bar(x + w/2, t_c, w, label='coupled (MaterialEngine)', color='#2980b9')

    # overhead annotation
    for i, r in enumerate(results):
        ax.annotate(f'{r["overhead"]:.1f}x',
                    xy=(x[i] + w/2, t_c[i]),
                    ha='center', va='bottom', fontsize=8, color='#2c3e50')

    ax.set_ylabel('wall time (s)')
    ax.set_title('Forward solve: static vs coupled')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha='right', fontsize=8)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'bench_forward_bar.png', dpi=150)
    plt.close(fig)


def plot_slsqp_bar(results):
    labels = [r['label'] for r in results]
    t_s = [r['t_static'] for r in results]
    t_c = [r['t_coupled'] for r in results]

    x = np.arange(len(labels))
    w = 0.35

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(x - w/2, t_s, w, label='static', color='#95a5a6')
    ax.bar(x + w/2, t_c, w, label='coupled', color='#c0392b')

    for i, r in enumerate(results):
        ax.annotate(f'{r["overhead"]:.1f}x',
                    xy=(x[i] + w/2, t_c[i]),
                    ha='center', va='bottom', fontsize=9, color='#2c3e50')

    ax.set_ylabel('wall time (s)')
    ax.set_title(f'SLSQP optimizer: static vs coupled ({results[0]["iters"]} iterations)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha='right', fontsize=9)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'bench_slsqp_bar.png', dpi=150)
    plt.close(fig)


def plot_scaling(fwd_results):
    """Overhead vs node count for the PICA 1L cases at different npl."""
    pica_runs = [r for r in fwd_results if r['label'].startswith('PICA 1L')]
    if len(pica_runs) < 2:
        return

    nodes = [r['nodes'] for r in pica_runs]
    overheads = [r['overhead'] for r in pica_runs]
    t_s = [r['t_static'] for r in pica_runs]
    t_c = [r['t_coupled'] for r in pica_runs]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Left: absolute times
    ax1.plot(nodes, t_s, 'o--', color='#95a5a6', lw=2, ms=7, label='static')
    ax1.plot(nodes, t_c, 's-',  color='#2980b9', lw=2, ms=7, label='coupled')
    ax1.set_xlabel('number of nodes')
    ax1.set_ylabel('wall time (s)')
    ax1.set_title('PICA 1-layer: wall time vs grid size')
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Right: overhead ratio
    ax2.plot(nodes, overheads, 'D-', color='#e74c3c', lw=2, ms=7)
    ax2.set_xlabel('number of nodes')
    ax2.set_ylabel('overhead (coupled / static)')
    ax2.set_title('PICA 1-layer: coupling overhead vs grid size')
    ax2.axhline(1.0, color='gray', ls='--', alpha=0.5)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / 'bench_scaling.png', dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    main()
