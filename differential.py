"""
Created on 03 17 2026

@author: RickFu
"""
import postprocessing as pp
import heatConduction as hc
import pandas as pd
import numpy as np
import os
import parameter as parameter

def main(para, cache, verbose=True):
    """ Adjoint calculation for spatially varying properties.

    Computes per-node gradients: dJ/dk, dJ/drho, dJ/dcp.
    The loss is J = 0.5 * (max_t T_L(t) - target)^2.

    Derivation:
        dJ/dp_j = -sum_{n=1}^{N} (w^n)^T * (dM/dp_j) * T^n
        The adjoint source is injected at t* = argmax_t T_L(t).

        dM/dk_j affects rows j-1, j, j+1 (interface conductivities)
        dM/drho_j and dM/dcp_j affect only row j:
            (dM/drho_j) = -(1/rho_j) * (M[j,:] - I[j,:])
            (dM/dcp_j)  = -(1/cp_j)  * (M[j,:] - I[j,:])

    Returns: dict with 'grad_k', 'grad_rho', 'grad_cp', 'lambda_profile', 'loss'
    """
    if verbose:
        print("Start Adjoint Calculation")
    target_temperature = para['back_wall_temperature_target']
    log = cache['Log']
    n_grid = cache['TProfile'][:, -1].size
    if verbose:
        print('Size of grid: ', n_grid)
    timeSteps = log.index
    if verbose:
        print('Time steps: ', len(timeSteps))

    # Variables
    lambda_current = np.zeros(n_grid)
    lambda_profile = np.zeros((len(timeSteps), n_grid))
    grad_k = np.zeros(n_grid)
    grad_rho = np.zeros(n_grid)
    grad_cp = np.zeros(n_grid)

    # Loss based on max backwall T over all timesteps
    T_backwall = cache['TProfile'][-1, :]
    t_star = np.argmax(T_backwall)
    T_max_bw = T_backwall[t_star]
    loss = 0.5 * (T_max_bw - target_temperature)**2
    # Adjoint source injected at t_star inside the backward loop (not terminal)

    # Properties (per-node arrays)
    rho = para['density']
    hcp = para['heatCapacity']
    dt = para['deltaTime']

    # BC info
    typeX0 = para['x=0 type']
    typeXL = para['x=L type']
    reradiate = para.get('re-radiate', False)
    if reradiate:
        sigma = para['stefanBoltzmann']
        eps_r = para['emissivity']
        T_amb = para['ambientTemperature']

    # Jacobian and per-node coefficients from forward solve
    # M_base is the Jacobian without re-radiation (constant across timesteps)
    M_base = cache['Jacobian'].copy()
    ce = cache['ce_arr']
    cw = cache['cw_arr']

    # Remove re-radiation from M_base so we can re-add it per timestep
    if reradiate:
        dx_arr = para['dx_array']
        T_last = cache['TProfile'][:, -1]
        if typeX0 == 'heatFlux':
            h_0 = 0.5 * (dx_arr[0] + dx_arr[0])
            rad0 = dt / (rho[0]*hcp[0]) * 2.0/h_0 * eps_r*sigma*4*T_last[0]**3
            M_base[0, 0] -= rad0
            M_base[0, 1] += rad0

    # Reverse time loop
    reversedtimeSteps = timeSteps[::-1]
    if verbose:
        print(' [Step]  [T_L]      [lambda_L]   [|grad_k|]')
    for ts in reversedtimeSteps:
        T_n = cache['TProfile'][:, ts]

        # Inject adjoint source at the timestep where backwall T is max
        if ts == t_star:
            lambda_current[-1] += (T_max_bw - target_temperature)

        # Rebuild M at this timestep (add re-radiation evaluated at T_n)
        if reradiate:
            M = M_base.copy()
            if typeX0 == 'heatFlux':
                h_0 = 0.5 * (dx_arr[0] + dx_arr[0])
                rad0 = dt / (rho[0]*hcp[0]) * 2.0/h_0 * eps_r*sigma*4*T_n[0]**3
                M[0, 0] += rad0
                M[0, 1] -= rad0
            MT = M.T
        else:
            MT = M_base.T

        # Backward propagation: solve M^T * w = lambda
        lambda_current = np.linalg.solve(MT, lambda_current)
        lam = lambda_current
        lambda_profile[ts-1, :] = lambda_current

        # --- grad_k: dM/dk_j affects rows j-1, j, j+1 ---
        # k_j enters k_{j-1/2} (east of row j-1, west of row j) with weight 1/2
        #       and k_{j+1/2} (east of row j, west of row j+1) with weight 1/2
        # At boundaries, ghost interface k = k[0] or k[-1] with weight 1.
        #
        # For interior j: (dM/dk_j @ T)_r contributions:
        #   row j-1: ce[j-1]/2 * (T[j-1] - T[j])
        #   row j:   cw[j]/2*(T[j]-T[j-1]) + ce[j]/2*(T[j]-T[j+1])
        #          = (cw[j]/2+ce[j]/2)*T[j] - cw[j]/2*T[j-1] - ce[j]/2*T[j+1]
        #   row j+1: cw[j+1]/2 * (T[j+1] - T[j])
        # Interior nodes (vectorized)
        grad_k[1:-1] -= (
            lam[:-2]  * ce[:-2]  / 2.0 * (T_n[:-2] - T_n[1:-1])
          + lam[1:-1] * (cw[1:-1] + ce[1:-1]) / 2.0 * T_n[1:-1]
          - lam[1:-1] * cw[1:-1] / 2.0 * T_n[:-2]
          - lam[1:-1] * ce[1:-1] / 2.0 * T_n[2:]
          + lam[2:]   * cw[2:]   / 2.0 * (T_n[2:] - T_n[1:-1])
        )
        # Boundary j=0: ghost interface dk/dk_0 = 1, east interface dk/dk_0 = 1/2
        #   row 0: cw[0]*1*T[0] + ce[0]/2*T[0] - (cw[0]*1 + ce[0]/2)*T_neighbor
        #   row 1: ce contribution from k_{1/2}: cw[1]/2*(T[1]-T[0])
        if typeX0 == 'heatFlux':
            v0_r0 = (cw[0] + ce[0]/2)*T_n[0] - (cw[0] + ce[0]/2)*T_n[1]
        elif typeX0 == 'fixedTemperature':
            v0_r0 = (cw[0] + ce[0]/2)*T_n[0] + (-ce[0]/2 + cw[0])*T_n[1]
        v0_r1 = cw[1]/2 * T_n[1] - cw[1]/2 * T_n[0]
        grad_k[0] -= lam[0] * v0_r0 + lam[1] * v0_r1
        # Boundary j=N-1: ghost interface dk/dk_{N-1} = 1, west dk = 1/2
        if typeXL == 'heatFlux':
            vN_rN = (ce[-1] + cw[-1]/2)*T_n[-1] - (ce[-1] + cw[-1]/2)*T_n[-2]
        elif typeXL == 'fixedTemperature':
            vN_rN = (ce[-1] + cw[-1]/2)*T_n[-1] + (ce[-1] - cw[-1]/2)*T_n[-2]
        vN_rN2 = ce[-2]/2 * T_n[-2] - ce[-2]/2 * T_n[-1]
        grad_k[-1] -= lam[-1] * vN_rN + lam[-2] * vN_rN2

        # --- grad_rho, grad_cp ---
        # At convergence: T^n - T^{n-1} = dt/(rho*cp) * D(T^n)
        # So dF_j/drho_j = (1/rho_j) * (T^n_j - T^{n-1}_j)
        #    dF_j/dcp_j  = (1/cp_j)  * (T^n_j - T^{n-1}_j)
        # This correctly handles boundary source terms.
        T_prev = cache['TProfile'][:, ts - 1]
        dT = T_n - T_prev
        grad_rho -= lam * dT / rho
        grad_cp  -= lam * dT / hcp

        if verbose:
            print(' [','{:3.0f}'.format(ts), ']',
                  ' [','{:8.2f}'.format(T_n[-1]),']',
                  ' [','{:10.4E}'.format(lambda_current[-1]),']',
                  ' [','{:10.4E}'.format(np.linalg.norm(grad_k)),']')

    if verbose:
        print('\nFinal loss:     ', '{:.6E}'.format(loss))
        print('|grad_k|:       ', '{:.6E}'.format(np.linalg.norm(grad_k)))
        print('|grad_rho|:     ', '{:.6E}'.format(np.linalg.norm(grad_rho)))
        print('|grad_cp|:      ', '{:.6E}'.format(np.linalg.norm(grad_cp)))

    return {
        'grad_k': grad_k, 'grad_rho': grad_rho, 'grad_cp': grad_cp,
        'lambda_profile': lambda_profile, 'loss': loss
    }


if __name__ == "__main__":
    para = parameter.main()
    outputDir = para['output']
    print('Output directory: ' + outputDir)
    os.makedirs(outputDir, exist_ok=True)
    results, cache = hc.solve(para)
    result = main(para, cache)