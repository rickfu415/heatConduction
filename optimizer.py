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
    k_layers = np.atleast_1d(np.asarray(para['layerConductivities'], dtype=float))
    rho_layers = np.atleast_1d(np.asarray(para['layerDensities'], dtype=float))
    cp_layers = np.atleast_1d(np.asarray(para['layerHeatCapacities'], dtype=float))
    t_layers = np.atleast_1d(np.asarray(para['layerThicknesses'], dtype=float))
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


def _rebuild_para(para_base, t):
    """Return a fresh para copy consistent with thickness vector t."""
    para = para_base.copy()
    t = np.asarray(t, dtype=float)
    para['layerThicknesses'] = t.copy()
    para['length'] = float(t.sum())
    para['material function'] = 'layered'
    parameter.normalize_conductivity(para)
    return para


def _layer_node_ranges(para):
    """Return list of (start_node, end_node) inclusive for each layer."""
    npl = para['nodesPerLayer']
    # np.atleast_1d guards against pandas coercing a 1-element array to a 0-d scalar
    n_layers = len(np.atleast_1d(para['layerConductivities']))
    return [(l * (npl - 1), l * (npl - 1) + (npl - 1)) for l in range(n_layers)]


def compute_thermal_sensitivities(para, eps_frac=0.002):
    """
    Run one forward solve + n_layers perturbations to compute thickness sensitivities.

    The adjoint (chain_rule) computes zero-sum boundary-shift derivatives that are
    correct for the redistribution optimizer (optimize_layered) but wrong for
    independent thickness variations needed by SLSQP. The correct independent
    dT/dt_l requires differentiating through ghost cells, BC terms, and non-uniform
    dx at layer interfaces — essentially the full forward assembly. Finite differences
    on the forward solver give the exact independent sensitivities at the cost of
    n_layers extra forward solves per call, which is negligible for 2–5 layers.

    Returns dict:
        'TProfile'      : np.ndarray (n_nodes, n_steps+1)
        'cache'         : forward solve cache
        'T_bw_max'      : float   — peak backwall temperature
        'T_layer_max'   : np.ndarray (n_layers,)  — peak T in each layer
        'dT_bw_dt'      : np.ndarray (n_layers,)  — dT_backwall_max / dt_l  (FD)
        'dT_layer_dt'   : np.ndarray (n_layers, n_layers)  — [i,j] = dT_layer_i_max / dt_j (FD)
    """
    # np.atleast_1d guards against pandas coercing a 1-element array to a 0-d scalar
    n_layers = len(np.atleast_1d(para['layerConductivities']))
    t_layers = np.atleast_1d(np.asarray(para['layerThicknesses'], dtype=float))

    # Base forward solve
    TProfile, fwd_cache = hc.solve(para, verbose=False)
    T_bw_max = float(np.max(TProfile[-1, :]))

    # Per-layer: find the node that attains the highest temperature over all time
    layer_ranges = _layer_node_ranges(para)
    T_layer_max = np.zeros(n_layers)
    obs_nodes = np.zeros(n_layers, dtype=int)
    for i, (s, e) in enumerate(layer_ranges):
        sub = TProfile[s:e+1, :]
        flat_idx = np.argmax(sub)
        node_local = np.unravel_index(flat_idx, sub.shape)[0]
        obs_nodes[i] = s + node_local
        T_layer_max[i] = float(TProfile[obs_nodes[i], :].max())

    # Thickness Jacobians via forward finite differences (one extra solve per layer).
    # Each perturbed solve uses the same obs_nodes as the base solve — the hottest
    # node location is assumed stable under small thickness perturbations.
    dT_bw_dt = np.zeros(n_layers)
    dT_layer_dt = np.zeros((n_layers, n_layers))
    for j in range(n_layers):
        eps = max(1e-5, eps_frac * t_layers[j])
        t_pert = t_layers.copy()
        t_pert[j] += eps
        TProfile_p, _ = hc.solve(_rebuild_para(para, t_pert), verbose=False)
        T_bw_p = float(np.max(TProfile_p[-1, :]))
        dT_bw_dt[j] = (T_bw_p - T_bw_max) / eps
        for i in range(n_layers):
            T_layer_p = float(TProfile_p[obs_nodes[i], :].max())
            dT_layer_dt[i, j] = (T_layer_p - T_layer_max[i]) / eps

    return {
        'TProfile': TProfile,
        'cache': fwd_cache,
        'T_bw_max': T_bw_max,
        'T_layer_max': T_layer_max,
        'dT_bw_dt': dT_bw_dt,
        'dT_layer_dt': dT_layer_dt,
    }


def optimize_mass_slsqp(
    para_base,
    t0,
    T_bw_limit,
    layer_service_temps,
    t_min,
    t_max=None,
    tol=1e-6,
    max_iter=300,
    callback=None,
):
    """ Constrained areal-mass minimization via SLSQP.

    Solves:
        minimize    m(t) = sum_i rho_i * t_i
        subject to  T_backwall_max(t) <= T_bw_limit
                    T_layer_i_max(t)  <= layer_service_temps[i]   (if finite)
                    t_i >= t_min[i]

    Gradients:
        - Mass objective gradient: analytic (rho_i).
        - Constraint Jacobians: adjoint-based via compute_thermal_sensitivities().

    Parameters
    ----------
    para_base : pd.Series
        Base parameter set. Never mutated; copies are made internally.
    t0 : np.ndarray, shape (n_layers,)
        Initial layer thicknesses (m). Should satisfy all constraints for reliable convergence.
    T_bw_limit : float
        Maximum allowed back wall temperature (K).
    layer_service_temps : list[float], length n_layers
        Per-layer service temperature limits (K). Pass np.inf to skip a layer.
    t_min : float or np.ndarray
        Lower bound(s) on layer thicknesses (m).
    t_max : float or np.ndarray or None
        Upper bound(s) on layer thicknesses (m). None means no upper bound.
    tol : float
        SLSQP optimality tolerance (ftol). Default 1e-6 matches scipy's SLSQP default.
    max_iter : int
        Maximum SLSQP iterations.
    callback : callable or None
        Called each iteration with a dict containing:
            iter, mass_kg_m2, t_layers_mm, T_bw_max, T_layer_max, constraint_violations

    Returns
    -------
    scipy.optimize.OptimizeResult
        Standard result plus extra attributes:
            .history   : list of per-iteration dicts
            .para_opt  : updated para Series at the optimum
    """
    from scipy.optimize import minimize, Bounds

    t0 = np.asarray(t0, dtype=float)
    n_layers = len(t0)

    rho = np.atleast_1d(np.asarray(para_base['layerDensities'], dtype=float))
    layer_service_temps = list(layer_service_temps)
    t_min_arr = np.broadcast_to(np.asarray(t_min, dtype=float), (n_layers,)).copy()

    # Auto-initialize: if t0 is deeply feasible (T_bw << T_bw_limit), scale t down
    # until T_bw ≈ T_bw_limit.  Starting near the constraint boundary avoids the
    # large initial QP step that SLSQP would take from a deeply-feasible point,
    # which can cause it to overshoot into the infeasible region and then climb back
    # to a wrong (heavier) local solution.
    _, _c0 = hc.solve(_rebuild_para(para_base, t0), verbose=False)
    _T_bw0 = float(np.max(_c0['TProfile'][-1, :]))
    if _T_bw0 < T_bw_limit - 10.0:
        _lo, _hi = 0.0, 1.0
        for _ in range(30):
            _mid = 0.5 * (_lo + _hi)
            _t_mid = np.maximum(_mid * t0, t_min_arr)
            _, _c_mid = hc.solve(_rebuild_para(para_base, _t_mid), verbose=False)
            _T_mid = float(np.max(_c_mid['TProfile'][-1, :]))
            if _T_mid < T_bw_limit:
                _hi = _mid
            else:
                _lo = _mid
        # _hi is smallest scale where T_bw < T_bw_limit; use it (feasible, near boundary)
        t0 = np.maximum(_hi * t0, t_min_arr)
        _, _c_new = hc.solve(_rebuild_para(para_base, t0), verbose=False)
        _T_new = float(np.max(_c_new['TProfile'][-1, :]))
        print(' Auto-init: t0 → {} mm  (T_bw: {:.1f}→{:.1f} K)'.format(
              (t0*1000).round(2).tolist(), _T_bw0, _T_new))

    # Variable scaling: work with s = t / t_scale so that x0 = 1 and BFGS steps are O(1).
    # Without scaling, BFGS (initialized to identity) produces steps of O(rho) in meter units,
    # which overshoot wildly on the first iteration for this nonlinear problem.
    t_scale = t0.copy()

    def _s2t(s):
        """Scaled → physical thicknesses."""
        return np.asarray(s, dtype=float) * t_scale

    # Closure-local eval cache: reuse forward solve when SLSQP calls f and g at same s
    _ecache = {'s_key': None, 'result': None}

    def _eval(s):
        s = np.asarray(s, dtype=float)
        if _ecache['s_key'] is None or not np.allclose(_ecache['s_key'], s, rtol=0, atol=1e-14):
            para = _rebuild_para(para_base, _s2t(s))
            _ecache['result'] = compute_thermal_sensitivities(para)
            _ecache['s_key'] = s.copy()
        return _ecache['result']

    # Objective in scaled space (gradient = rho * t_scale)
    rho_scaled = rho * t_scale

    def f(s):
        return float(np.dot(rho_scaled, s))

    def jac_f(s):
        return rho_scaled.copy()

    # Constraints in scaled space. Jacobians pick up the t_scale factor via chain rule:
    #   dg/ds_j = dg/dt_j * t_scale_j
    # Also scaled by reference T for additional SLSQP conditioning.
    constraints = []
    constraints.append({
        'type': 'ineq',
        'fun': lambda s: (T_bw_limit - _eval(s)['T_bw_max']) / T_bw_limit,
        'jac': lambda s: -_eval(s)['dT_bw_dt'] * t_scale / T_bw_limit,
    })
    for i in range(n_layers):
        T_svc = float(layer_service_temps[i])
        if not np.isinf(T_svc):
            constraints.append({
                'type': 'ineq',
                'fun': lambda s, _i=i, _T=T_svc: (_T - _eval(s)['T_layer_max'][_i]) / _T,
                'jac': lambda s, _i=i, _T=T_svc: -_eval(s)['dT_layer_dt'][_i, :] * t_scale / _T,
            })

    bounds_scaled = Bounds(
        lb=t_min_arr / t_scale,
        ub=np.full(n_layers, np.inf) if t_max is None
           else np.asarray(t_max, dtype=float) / t_scale,
        keep_feasible=True,
    )

    # History + callback (convert s back to t for the user)
    history = []
    iter_counter = [0]

    def slsqp_callback(s):
        iter_counter[0] += 1
        t = _s2t(s)
        r = _eval(s)
        entry = {
            'iter': iter_counter[0],
            'mass_kg_m2': float(np.dot(rho, t)),
            't_layers_mm': (t * 1000).tolist(),
            'T_bw_max': r['T_bw_max'],
            'T_layer_max': r['T_layer_max'].tolist(),
            'constraint_violations': {
                'bw': max(0.0, r['T_bw_max'] - T_bw_limit),
                'layer': [
                    max(0.0, r['T_layer_max'][i] - layer_service_temps[i])
                    for i in range(n_layers)
                    if not np.isinf(layer_service_temps[i])
                ],
            },
        }
        history.append(entry)
        if callback is not None:
            callback(entry)

    print('='*70)
    print(' Mass Minimization: {} layers  |  T_bw_limit={:.0f} K'.format(n_layers, T_bw_limit))
    print(' Initial mass: {:.3f} kg/m²  |  t0: {} mm'.format(
          float(np.dot(rho, t0)), (t0 * 1000).round(3).tolist()))
    print('='*70)

    s0 = np.ones(n_layers)   # x0 in scaled space (t0 / t_scale = 1)

    opt_result = minimize(
        fun=f,
        x0=s0,
        jac=jac_f,
        method='SLSQP',
        bounds=bounds_scaled,
        constraints=constraints,
        options={'maxiter': max_iter, 'ftol': tol, 'disp': True},
        callback=slsqp_callback,
    )

    # Convert back to physical thicknesses
    t_opt = _s2t(opt_result.x)
    opt_result.x = t_opt
    opt_result.fun = float(np.dot(rho, t_opt))
    opt_result.history = history
    opt_result.para_opt = _rebuild_para(para_base, t_opt)

    print('='*70)
    print(' SLSQP status: {}'.format(opt_result.message))
    print(' Optimized mass: {:.3f} kg/m²  |  t_opt: {} mm'.format(
          float(opt_result.fun), (opt_result.x * 1000).round(3).tolist()))
    print('='*70)

    return opt_result


def optimize_layered(para, max_iter=200, tol=1.0, callback=None):
    """ Gradient descent optimizer for layered TPS.

    Optimizes layer thicknesses using adjoint gradients + chain rule.
    Material properties (k, rho, cp) per layer are fixed.
    Total thickness is free to change to meet the back wall T target.

    Return: optimized parameter series, history list
    """
    target = para['back_wall_temperature_target']
    # np.atleast_1d guards against pandas coercing 1-element arrays to 0-d scalars
    k_layers = np.atleast_1d(np.asarray(para['layerConductivities'], dtype=float)).copy()
    t_layers = np.atleast_1d(np.asarray(para['layerThicknesses'], dtype=float)).copy()
    lr = para.get('learning_rate', 1e-9)
    n_layers = len(k_layers)
    t_min = max(1e-5, 0.01 * t_layers.min())

    # Edge case: check if redistribution alone can meet target.
    # Best redistribution: maximize insulator thickness, minimize others at t_min.
    # Use per-layer 10% minimums so the hot-face layer never becomes critically thin
    # (a global t_min driven by a thin back-wall layer would make the hot-face
    # layer numerically unstable under high heat flux).
    best_insulator = np.argmin(k_layers)
    t_min_arr = np.maximum(1e-4, 0.1 * t_layers)
    t_check = t_min_arr.copy()
    t_check[best_insulator] = t_layers.sum() - t_min_arr.sum() + t_min_arr[best_insulator]
    para_check = para.copy()
    para_check['layerThicknesses'] = t_check
    para_check['length'] = float(t_check.sum())
    para_check['material function'] = 'layered'
    T_best = None
    try:
        _, cache_check = hc.solve(para_check, verbose=False)
        T_best = np.max(cache_check['TProfile'][-1, :])
        needs_growth = T_best > target
        print(' Best redistribution check: max T_L = {:.1f} K (layer {} at {:.3f}mm, others at t_min)'.format(
              T_best, best_insulator, t_check[best_insulator]*1000))
        if needs_growth:
            print(' -> Current total thickness insufficient, allowing growth')
        else:
            print(' -> Redistribution sufficient, keeping total thickness constant')
    except Exception:
        # Pre-check numerically unstable (e.g. thin hot-face layer under high heat flux).
        # Default to needs_growth=True — safe conservative choice.
        needs_growth = True
        T_best = None
        print(' Pre-check solve failed (thin layer / high flux); assuming growth needed.')

    # Emit pre-check result to UI before iteration starts
    if callback:
        callback({
            'type': 'precheck',
            'best_insulator_layer': int(best_insulator),
            'T_best': float(T_best) if T_best is not None else None,
            'needs_growth': bool(needs_growth),
            't_check_mm': (t_check * 1000).tolist(),
            'target': float(target),
        })

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

        # Forward solve — guard against numerical blow-up from extreme thin layers
        try:
            _, cache = hc.solve(para, verbose=False)
        except Exception as e:
            print(' Forward solve failed at iteration {} ({}). Stopping early.'.format(iteration, e))
            break
        T_L = np.max(cache['TProfile'][-1, :])  # max backwall T over time
        loss = 0.5 * (T_L - target)**2

        # Adjoint → per-node gradients for k, rho, cp
        adjoint_result = diff.main(para, cache, verbose=False)

        # Chain rule → thickness gradient (redistribution only, sums to ~0)
        grad_t = chain_rule(adjoint_result, para)
        # When growth is needed and T_L exceeds target, add a uniform growth
        # component proportional to each layer's insulating value (1/k).
        # This breaks the zero-sum symmetry so total thickness can increase.
        # Once T_L drops below target, disable growth permanently (sufficient thickness reached).
        if needs_growth and T_L <= target:
            needs_growth = False
            print(' -> Thickness now sufficient, switching to redistribution only')
        if needs_growth:
            growth_dir = (1.0 / k_layers)
            growth_dir = growth_dir / np.linalg.norm(growth_dir)
            grad_norm_now = np.linalg.norm(grad_t)
            if grad_norm_now == 0.0:
                # No boundary gradients (e.g. single-layer): grow at the clip
                # limit (50%/step) until T_L drops below target.
                grad_t = -growth_dir / lr
            else:
                grad_t -= growth_dir * (T_L - target) * grad_norm_now
        elif np.linalg.norm(grad_t) == 0.0:
            # Chain rule gave zero gradient (e.g. single layer, no boundaries).
            # Estimate dJ/dt via a 1% forward finite-difference perturbation.
            dt_fd = max(1e-5, 0.01 * float(t_layers.sum()))
            t_fd = np.atleast_1d(t_layers) + dt_fd
            para_fd = para.copy()
            para_fd['layerThicknesses'] = t_fd
            para_fd['length'] = float(t_fd.sum())
            try:
                _, cache_fd = hc.solve(para_fd, verbose=False)
                T_L_fd = float(np.max(cache_fd['TProfile'][-1, :]))
                dT_dt = (T_L_fd - T_L) / dt_fd  # dT_L / dt  (< 0: thicker → cooler)
                if dT_dt != 0.0:
                    dJ_dt = (T_L - target) * dT_dt
                    grad_t = np.full(n_layers, dJ_dt / n_layers)
            except Exception:
                pass  # keep grad_t = 0; iteration stalls but does not crash
        grad_norm = np.linalg.norm(grad_t)

        # Log
        history.append({
            'iter': iteration, 'T_L': T_L, 'loss': loss,
            'grad_norm': grad_norm, 't_layers': t_layers.copy()
        })
        t_str = np.array2string(t_layers*1000, precision=3, separator=',')
        print(' {:4d}  {:10.2f}  {:12.4E}  {:12.4E}  {:>30s}'.format(
              iteration, T_L, loss, grad_norm, t_str + ' mm'))

        if callback is not None:
            callback({
                'iter': iteration, 'T_L': float(T_L), 'loss': float(loss),
                'grad_norm': float(grad_norm),
                't_layers': t_layers.copy().tolist(),
                'total_thickness': float(t_layers.sum()),
            })

        if loss < tol:
            print('-'*70)
            print(' Converged! Loss {:.4E} < tol {:.4E}'.format(loss, tol))
            break

        # Gradient descent with fixed learning rate
        step = lr * grad_t
        # Limit step: no layer can change by more than 50% per iteration
        max_step = 0.5 * t_layers
        step = np.clip(step, -max_step, max_step)
        t_new = t_layers - step
        t_new = np.maximum(t_new, t_min)
        # If thickness is already sufficient and there are multiple layers,
        # only redistribute (keep total constant). For a single layer there is
        # nothing to redistribute — allow total thickness to change.
        if not needs_growth and n_layers > 1:
            t_new = t_new * t_layers.sum() / t_new.sum()
        t_layers = t_new
        para['length'] = float(t_layers.sum())

    print('='*70)
    print(' Fixed k_layers:', k_layers)
    print(' Optimized t_layers (mm):', t_layers*1000)
    T_L_final = history[-1]['T_L'] if history else float('nan')
    print(' Final T_L = {:.4f} K  (target = {:.1f} K)'.format(T_L_final, target))
    print('='*70)

    para['layerThicknesses'] = t_layers
    para['length'] = float(t_layers.sum())
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
