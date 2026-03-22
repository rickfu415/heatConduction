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
    """
    Main function

    """
    if verbose:
        print("Start Adjoint Calculation")
    target_temperature = para['back_wall_temperature_target']
    log = cache['Log']
    T_final_distribution = cache['TProfile'][:, -1]
    n_grid = T_final_distribution.size
    if verbose:
        print('Size of grid: ', n_grid)
    timeSteps = log.index
    if verbose:
        print('Time steps: ', len(timeSteps))

    # define variables
    lambda_current = np.zeros(n_grid)
    lambda_profile = np.zeros((len(timeSteps), n_grid))
    grad_k = 0.0

    # Terminal condition: dJ/dT^N (derivative of loss, not loss itself)
    T_final = cache['TProfile'][:, -1]
    loss = 0.5 * (T_final[-1] - target_temperature)**2
    lambda_current[-1] = T_final[-1] - target_temperature

    # Build the Jacobian M (constant for constant properties)
    # Reuse assemble to get M at the final state
    k = para['conductivity']
    rho = para['density']
    hcp = para['heatCapacity']
    dt = para['deltaTime']
    length = para['length']
    alpha = k / rho / hcp
    dx = length / (n_grid - 1)

    # dM/dk: derivative of Jacobian w.r.t. conductivity
    # M_ii = 1 + 2*alpha*dt/dx^2, M_{i,i+-1} = -alpha*dt/dx^2
    # alpha = k/(rho*cp), so d(alpha)/dk = 1/(rho*cp)
    dalpha_dk = 1.0 / (rho * hcp)
    dM_dk = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        dM_dk[i, i] = 2 * dalpha_dk * dt / dx**2
        if i > 0:
            dM_dk[i, i-1] = -dalpha_dk * dt / dx**2
        if i < n_grid - 1:
            dM_dk[i, i+1] = -dalpha_dk * dt / dx**2

    # Reconstruct M from the forward solve Jacobian
    M = cache['Jacobian']
    MT = M.T

    # Reverse time loop for adjoint calculation
    reversedtimeSteps = timeSteps[::-1]
    if verbose:
        print(' [Step]  [T_L]      [lambda_L]   [grad_k]')
    for ts in reversedtimeSteps:
        T_n = cache['TProfile'][:, ts]

        # Accumulate gradient: dJ/dk += lambda^T * (dM/dk * T^n)
        grad_k += lambda_current @ (dM_dk @ T_n)

        # Store adjoint
        lambda_profile[ts-1, :] = lambda_current

        # Backward propagation: solve M^T * lambda^{n-1} = lambda^n
        lambda_current = np.linalg.solve(MT, lambda_current)

        # output results
        if verbose:
            print(' [','{:3.0f}'.format(ts), ']',
                  ' [','{:8.2f}'.format(T_n[-1]),']',
                  ' [','{:10.4E}'.format(lambda_current[-1]),']',
                  ' [','{:10.4E}'.format(grad_k),']')

    if verbose:
        print('\nFinal loss:    ', '{:.6E}'.format(loss))
        print('Gradient dJ/dk:', '{:.6E}'.format(grad_k))
    return grad_k, lambda_profile


if __name__ == "__main__":
    para = parameter.main()
    outputDir = para['output']
    print('Output directory: ' + outputDir)
    os.makedirs(outputDir, exist_ok=True)
    results, cache = hc.solve(para)
    grad_k, lambda_profile = main(para, cache)