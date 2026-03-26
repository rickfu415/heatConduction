"""
Created on 03 21 2026

@author: RickFu

Adjoint-based optimizer for layered TPS heat conduction.

Framework:
    forward solver → adjoint (dJ/dk, dJ/drho, dJ/dcp per node)
    → chain rule (map to layer thickness gradient) → gradient descent
"""
import heatConduction as hc
import differential as diff
import postprocessing as pp
import parameter
import numpy as np
import os


def chain_rule(adjoint_result, para):
    """ Map per-node adjoint gradients to layer thickness gradient.

    When a layer boundary moves right by delta, nodes switch from
    layer l+1 to layer l. Their k, rho, cp all change. The total
    sensitivity is the sum of contributions from all three properties.

    dJ/dt_l = [dJ/dk   * (k_l   - k_{l+1})
             + dJ/drho * (rho_l - rho_{l+1})
             + dJ/dcp  * (cp_l  - cp_{l+1})] * grad_at_boundary / dx

    Returns: grad_t_layers (length n_layers)
    """
    k_layers = para['layerConductivities']
    rho_layers = para['layerDensities']
    cp_layers = para['layerHeatCapacities']
    t_layers = para['layerThicknesses']
    npl = para['nodesPerLayer']
    n_layers = len(k_layers)

    grad_k = adjoint_result['grad_k']
    grad_rho = adjoint_result['grad_rho']
    grad_cp = adjoint_result['grad_cp']
    dx_arr = para['dx_array']

    grad_t = np.zeros(n_layers)
    for l in range(n_layers - 1):
        # With layered grid, boundary between layer l and l+1 is at a known node
        bnd_node = (l + 1) * (npl - 1)
        # dx at the boundary (average of both sides)
        dx_local = 0.5 * (dx_arr[bnd_node - 1] + dx_arr[bnd_node])

        sensitivity = (
            (k_layers[l]   - k_layers[l+1])   * grad_k[bnd_node]
          + (rho_layers[l] - rho_layers[l+1]) * grad_rho[bnd_node]
          + (cp_layers[l]  - cp_layers[l+1])  * grad_cp[bnd_node]
        ) / dx_local

        grad_t[l] += sensitivity
        grad_t[l + 1] -= sensitivity

    return grad_t


def optimize_layered(para, max_iter=200, tol=1.0):
    """ Gradient descent optimizer for layered TPS.

    Optimizes layer thicknesses using adjoint gradients + chain rule.
    Material properties (k, rho, cp) per layer are fixed.

    Return: optimized parameter series, history list
    """
    target = para['back_wall_temperature_target']
    k_layers = para['layerConductivities'].copy()
    t_layers = para['layerThicknesses'].copy()
    L = para['length']
    n_layers = len(k_layers)
    t_min = L * 0.01

    history = []

    print('='*70)
    print(' Layered TPS Optimization: {} layers'.format(n_layers))
    print(' Target T_L = {:.1f} K'.format(target))
    print(' Fixed k_layers:', k_layers)
    print(' Initial t_layers:', t_layers)
    print('='*70)
    print(' {:>4s}  {:>10s}  {:>12s}  {:>12s}  {:>30s}'.format(
          'Iter', 'T_L', 'Loss', '|grad_t|', 't_layers'))
    print('-'*70)

    for iteration in range(1, max_iter + 1):
        para['layerThicknesses'] = t_layers
        para['material function'] = 'layered'

        # Forward solve
        _, cache = hc.solve(para, verbose=False)
        T_L = np.max(cache['TProfile'][-1, :])  # max backwall T over time
        loss = 0.5 * (T_L - target)**2

        # Adjoint → per-node gradients for k, rho, cp
        adjoint_result = diff.main(para, cache, verbose=False)

        # Chain rule → thickness gradient
        grad_t = chain_rule(adjoint_result, para)
        grad_norm = np.linalg.norm(grad_t)

        # Log
        history.append({
            'iter': iteration, 'T_L': T_L, 'loss': loss,
            'grad_norm': grad_norm, 't_layers': t_layers.copy()
        })
        t_str = np.array2string(t_layers*1000, precision=3, separator=',')
        print(' {:4d}  {:10.2f}  {:12.4E}  {:12.4E}  {:>30s}'.format(
              iteration, T_L, loss, grad_norm, t_str + ' mm'))

        if loss < tol:
            print('-'*70)
            print(' Converged! Loss {:.4E} < tol {:.4E}'.format(loss, tol))
            break

        # Gradient descent: step size adapts with loss
        # Large steps early (5% of L), shrinking as loss decreases
        max_step = 0.05 * L * min(1.0, loss / 1000.0)
        max_step = max(max_step, 0.001 * L)  # floor at 0.1%
        step = grad_t / (grad_norm + 1e-30) * max_step
        t_layers = t_layers - step
        t_layers = np.maximum(t_layers, t_min)
        t_layers = t_layers * L / t_layers.sum()

    print('='*70)
    print(' Fixed k_layers:', k_layers)
    print(' Optimized t_layers (mm):', t_layers*1000)
    print(' Final T_L = {:.4f} K  (target = {:.1f} K)'.format(T_L, target))
    print('='*70)

    para['layerThicknesses'] = t_layers
    return para, history


if __name__ == "__main__":
    para = parameter.main()
    outputDir = para['output']
    os.makedirs(outputDir, exist_ok=True)

    para, history = optimize_layered(para)

    # Final forward solve for postprocessing
    results, cache = hc.solve(para)
    T = pp.preprocess(para, results)
    T.to_csv(os.path.join(outputDir, 'solutionHistory_optimized.csv'))

    import pandas as pd
    df_hist = pd.DataFrame(history)
    df_hist.to_csv(os.path.join(outputDir, 'optimizationHistory.csv'), index=False)

    pp.evolutionField(T, outputDir)
    positions = pp.probePositions(para)
    pp.thermalCouplePlot(T, positions, outputDir)
    print('\nResults saved to:', outputDir)
