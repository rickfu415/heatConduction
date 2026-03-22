"""
Created on 03 21 2026

@author: RickFu

Gradient-based optimizer for heat conduction inverse problem.
Finds the conductivity k such that the back-wall temperature
matches the target at the final time.
"""
import heatConduction as hc
import differential as diff
import postprocessing as pp
import parameter
import os


def optimize(para, max_iter=200, tol=1e-2):
    """ Newton's method optimization loop

    Iterates forward solve -> adjoint -> Newton update on k until
    the back-wall temperature matches the target.

    Uses: dT_L/dk = dJ/dk / (T_L - target)
    Newton step: k_new = k - (T_L - target) / (dT_L/dk)
                       = k - (T_L - target)^2 / (dJ/dk)

    Inputs:
        para: parameter series
        max_iter: maximum optimization iterations
        tol: convergence tolerance on loss

    Return: optimized parameter series, history list
    """

    target = para['back_wall_temperature_target']
    k = para['conductivity']
    history = []
    k_prev = None
    T_L_prev = None

    print('='*70)
    print(' Optimization: find k so that T_L(t_final) = {:.1f} K'.format(target))
    print(' Initial k = {:.4f} W/(m*K)'.format(k))
    print('='*70)
    print(' {:>4s}  {:>12s}  {:>10s}  {:>12s}  {:>12s}'.format(
          'Iter', 'k', 'T_L', 'Loss', 'dT_L/dk'))
    print('-'*70)

    for iteration in range(1, max_iter + 1):
        # Update conductivity in parameters
        para['conductivity'] = k

        # Forward solve (silent)
        _, cache = hc.solve(para, verbose=False)

        # Back-wall temperature at final time
        T_L = cache['TProfile'][-1, -1]
        loss = 0.5 * (T_L - target)**2

        # Estimate dT_L/dk:
        #   Iter 1: use adjoint gradient (local linearization)
        #   Iter 2+: use secant method (actual finite difference)
        if k_prev is not None and abs(k - k_prev) > 1e-12:
            dTL_dk = (T_L - T_L_prev) / (k - k_prev)
        else:
            grad_k, _ = diff.main(para, cache, verbose=False)
            dTL_dk = grad_k / (T_L - target)

        # Log
        history.append({'iter': iteration, 'k': k,
                        'T_L': T_L, 'loss': loss, 'dTL_dk': dTL_dk})
        print(' {:4d}  {:12.4f}  {:10.2f}  {:12.4E}  {:12.4E}'.format(
              iteration, k, T_L, loss, dTL_dk))

        # Check convergence
        if loss < tol:
            print('-'*70)
            print(' Converged! Loss {:.4E} < tol {:.4E}'.format(loss, tol))
            break

        # Store previous values for secant
        k_prev = k
        T_L_prev = T_L

        # Newton/secant step: k_new = k - (T_L - target) / (dT_L/dk)
        k = k - (T_L - target) / dTL_dk

        # Safety: keep k positive
        if k <= 0:
            k = 1e-3
            print(' Warning: k clamped to {:.4f}'.format(k))

    print('='*70)
    print(' Optimized k = {:.6f} W/(m*K)'.format(k))
    print(' Final T_L   = {:.4f} K  (target = {:.1f} K)'.format(T_L, target))
    print('='*70)

    para['conductivity'] = k
    return para, history


if __name__ == "__main__":
    para = parameter.main()
    outputDir = para['output']
    os.makedirs(outputDir, exist_ok=True)

    # Run optimizer
    para, history = optimize(para)

    # Final forward solve with optimized k for postprocessing
    results, cache = hc.solve(para)
    T = pp.preprocess(para, results)
    T.to_csv(os.path.join(outputDir, 'solutionHistory_optimized.csv'))

    # Save optimization history
    import pandas as pd
    df_hist = pd.DataFrame(history)
    df_hist.to_csv(os.path.join(outputDir, 'optimizationHistory.csv'), index=False)

    # Save figures
    pp.evolutionField(T, outputDir)
    positions = [0, 0.002, 0.004, 0.006, 0.008, 0.01]
    pp.thermalCouplePlot(T, positions, outputDir)
    times = [0, 2, 4, 6, 8, 10]
    pp.temperatureDistribution(T, times, outputDir)

    print('\nResults saved to:', outputDir)
