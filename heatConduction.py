# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:15:28 2019

@author: RickFu
"""
import numpy as np
import pandas as pd
from scipy.linalg import solve_banded
import parameter
import utility
import time


def assemble(para, cache):
    """ Assemble linear system Jacobian * dx = F

    Supports spatially varying conductivity k(x) as a per-node array.
    Uses the divergence form: d/dx[k(x) * dT/dx] with arithmetic-mean
    interface conductivities k_{i+1/2} = (k[i] + k[i+1]) / 2.

    Return: dictionary containing cache data
    """

    k = para['conductivity']      # numpy array of length numberOfNode
    rho = para['density']          # numpy array of length numberOfNode
    hcp = para['heatCapacity']     # numpy array of length numberOfNode
    dt = para['deltaTime']
    numberOfNode = para['numberOfNode']
    dx_arr = para['dx_array']      # non-uniform grid spacing (length N-1)

    # BC informations
    typeX0 = para['x=0 type']
    valueX0 = para['x=0 value']
    typeXL = para['x=L type']
    valueXL = para['x=L value']
    reradiate = para.get('re-radiate', False)

    # Re-radiation parameters
    if reradiate:
        sigma = para['stefanBoltzmann']
        eps_r = para['emissivity']
        T_amb = para['ambientTemperature']

    # Containers
    T = cache['T']; T0 = cache['T0']
    F = cache['F']; Jacobian = cache['Jacobian']

    # Precompute interface conductivities and spacings (vectorized)
    N = numberOfNode
    k_half = np.empty(N + 1)
    dx_half = np.empty(N + 1)
    k_half[1:N] = 0.5 * (k[:-1] + k[1:])
    dx_half[1:N] = dx_arr
    k_half[0] = k[0];    dx_half[0] = dx_arr[0]
    k_half[N] = k[-1];   dx_half[N] = dx_arr[-1]

    # Compute effective heat flux at boundaries (including re-radiation)
    qX0 = valueX0
    qXL = valueXL
    if reradiate:
        if typeX0 == 'heatFlux':
            qX0 = valueX0 - eps_r * sigma * (T[0, 0]**4 - T_amb**4)

    # Boundary ghost values for temperature
    if typeX0 == 'heatFlux':
        Ug1 = utility.fixedGradient(qX0, k[0], dx_arr[0], T[1])
    elif typeX0 == 'fixedTemperature':
        Ug1 = utility.fixedValue(valueX0, T[1])

    if typeXL == 'heatFlux':
        Ug2 = utility.fixedGradient(qXL, k[-1], dx_arr[-1], T[-2])
    elif typeXL == 'fixedTemperature':
        Ug2 = utility.fixedValue(valueXL, T[-2])

    # Vectorized Jacobian assembly
    ke = k_half[1:]   # east interface conductivities (length N)
    kw = k_half[:-1]  # west interface conductivities (length N)
    dxe = dx_half[1:]
    dxw = dx_half[:-1]
    h = 0.5 * (dxw + dxe)

    ce_arr = dt / (rho * hcp * h * dxe)
    cw_arr = dt / (rho * hcp * h * dxw)

    # Diagonal: 1 + ce*k_east + cw*k_west
    diag = 1.0 + ce_arr * ke + cw_arr * kw

    # Off-diagonals (interior nodes)
    upper = -ce_arr[:-1] * ke[:-1]  # Jacobian[i][i+1] for i=0..N-2
    lower = -cw_arr[1:] * kw[1:]    # Jacobian[i][i-1] for i=1..N-1

    # Boundary corrections for off-diagonals
    if typeX0 == 'heatFlux':
        upper[0] = -(ce_arr[0] * ke[0] + cw_arr[0] * kw[0])
    elif typeX0 == 'fixedTemperature':
        upper[0] = -ce_arr[0] * ke[0] + cw_arr[0] * kw[0]

    if typeXL == 'heatFlux':
        lower[-1] = -(ce_arr[-1] * ke[-1] + cw_arr[-1] * kw[-1])
    elif typeXL == 'fixedTemperature':
        lower[-1] = ce_arr[-1] * ke[-1] - cw_arr[-1] * kw[-1]

    # Re-radiation Jacobian correction at x=0
    if reradiate:
        if typeX0 == 'heatFlux':
            h_0 = 0.5 * (dx_half[0] + dx_half[1])
            rad_jac_0 = dt / (rho[0]*hcp[0]) * 2.0 / h_0 * eps_r * sigma * 4 * T[0, 0]**3
            diag[0] += rad_jac_0
            upper[0] -= rad_jac_0

    # Fill Jacobian matrix
    Jacobian[:] = 0
    np.fill_diagonal(Jacobian, diag)
    np.fill_diagonal(Jacobian[:-1, 1:], upper)
    np.fill_diagonal(Jacobian[1:, :-1], lower)

    # Calculate F using variable-coefficient diffusion with non-uniform grid
    diffusion = utility.variableCoefficientDiffusion(T, k, dx_arr, Ug1, Ug2)
    rho_cp = (rho * hcp).reshape(-1, 1)
    F = T - T0 - dt / rho_cp * diffusion

    # Store in cache (ce/cw arrays used by adjoint for grad_k)
    cache['F'] = -F; cache['Jacobian'] = Jacobian
    cache['ce_arr'] = ce_arr; cache['cw_arr'] = cw_arr
    return cache


def initialize(para):
    """ Initialize key data
    
    T: current step temperature
    T0: last step temperature
    TProfile: temperature results in time and space
    F: B as right hand side of Ax = B
    Jacobian: A as left had side of Ax = B
    
    Return: a dictionary
    """
    
    numberOfNode = para['numberOfNode']
    numOfTimeStep = para['numberOfTimeStep']
    Tic = para['IC value']
    T = np.full((numberOfNode, 1), Tic)
    T0 = np.full((numberOfNode, 1), Tic)
    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    F = np.zeros((numberOfNode, 1))
    Jacobian = np.zeros((numberOfNode, numberOfNode))
    TProfile[:,0] = T.reshape(1,-1)
    cache = {'T':T,'T0':T0,'TProfile':TProfile,
             'F':F,'Jacobian':Jacobian,
             'Log':pd.DataFrame()}
    return cache


def solveLinearSystem(para, cache):
    """ Solve Ax=B using tridiagonal (banded) solver — O(N) instead of O(N³).

    The Jacobian is tridiagonal: pack into banded form for scipy.solve_banded.
    """
    relax = para['relaxation']
    A = cache['Jacobian']
    B = cache['F']
    N = A.shape[0]
    # Pack tridiagonal into banded form: ab[0] = upper, ab[1] = diag, ab[2] = lower
    ab = np.zeros((3, N))
    ab[0, 1:] = np.diag(A, 1)   # upper diagonal
    ab[1, :]  = np.diag(A, 0)   # main diagonal
    ab[2, :-1] = np.diag(A, -1) # lower diagonal
    dT = solve_banded((1, 1), ab, B.ravel()).reshape(-1, 1)
    T = cache['T']
    T = dT * relax + T
    cache['T'] = T
    cache['dT'] = dT
    return cache


def storeUpdateResult(cache):
    """ Store results
    Update T0
    Store temperaure results into a dataframe and 
    save it in the cache.
    """
    
    timeStep = cache['ts']
    TProfile = cache['TProfile']
    T = cache['T']
    cache['T0'] = T.copy()
    TProfile[:,timeStep] = T.reshape(1,-1)
    return cache


def newtonIteration(para, cache, verbose=True):
    """ Newton's Iteration for Equation System

    Process:
        1. Get max iteratino, convergence limit
        2. Call assemble function to get Jacobian and F(RHS)
        3. Solve for dT, update solution
        4. Evaluate F, get value of 2-norm
        5. If solution converged, break, output to screen and
           return cache.

    """

    maxIteration = para['maxIteration']
    convergence = para['convergence']
    dt = para['deltaTime']
    log = cache['Log']
    ts = cache['ts']
    for n in range(maxIteration):
        cache = assemble(para, cache)
        F = cache['F']
        norm = np.linalg.norm(F)
        if norm < convergence:
            log.loc[ts,'PhysicalTime'] = dt*ts
            log.loc[ts,'Iteration'] = n+1
            log.loc[ts,'Residual'] = norm
            break
        cache = solveLinearSystem(para, cache)
    print_freq = int(para.get('print_frequency', 1))
    if verbose and ts % print_freq == 0:
        T = cache['T']
        print(' [','{:3.0f}'.format(ts), ']',
              ' [','{:6.2f}'.format(ts*dt),']',
              ' [','{:2.0f}'.format(n+1), ']',
              ' [','{:8.2E}'.format(norm),']',
              ' [','{:8.2f}'.format(T[0,0]),']',
              ' [','{:8.2f}'.format(T[-1,0]),']')
    return cache


def solve(para, verbose=True):
    """ Main function to solve heat conduction

    Input: a Pandas series containing all parameters

    Process:
        1. Initialize cache
        2. Time marching
        3. Newton's iteration for discretized PDE for singe time
           step
        4. Update T, save result to T profile

    Return: temperature profile as final result
    """

    if verbose:
        print(" Heat Conduction Solver")
    start = time.time()
    para = parameter.normalize_conductivity(para)
    cache = initialize(para)
    numOfTimeStep = para['numberOfTimeStep']
    dt = para['deltaTime']

    # Detect time-varying BCs (2D array [[t0,q0],[t1,q1],...])
    val_x0 = para['x=0 value']
    val_xL = para['x=L value']
    flux_profile_x0 = None
    flux_profile_xL = None
    if isinstance(val_x0, np.ndarray) and val_x0.ndim == 2:
        flux_profile_x0 = val_x0.copy()
    if isinstance(val_xL, np.ndarray) and val_xL.ndim == 2:
        flux_profile_xL = val_xL.copy()

    if verbose:
        print(' [Step] [Pysical Time] [Iteration] [Residue]    [T_0]      [T_L]')
    flux_history_x0 = np.zeros(numOfTimeStep + 1)
    for timeStep in range(1, numOfTimeStep+1):
        cache['ts'] = timeStep
        t_phys = timeStep * dt
        if flux_profile_x0 is not None:
            para['x=0 value'] = float(np.interp(t_phys, flux_profile_x0[:, 0], flux_profile_x0[:, 1]))
        if flux_profile_xL is not None:
            para['x=L value'] = float(np.interp(t_phys, flux_profile_xL[:, 0], flux_profile_xL[:, 1]))
        flux_history_x0[timeStep] = para['x=0 value'] if np.isscalar(para['x=0 value']) else float(para['x=0 value'])
        cache = newtonIteration(para, cache, verbose=verbose)
        cache = storeUpdateResult(cache)

    # Store flux history for potential adjoint use; restore original BC values
    cache['flux_history_x0'] = flux_history_x0
    if flux_profile_x0 is not None:
        para['x=0 value'] = flux_profile_x0
    if flux_profile_xL is not None:
        para['x=L value'] = flux_profile_xL
    TProfile = cache['TProfile']
    runtime = time.time() - start
    if verbose:
        print('[Cost] CPU time spent','%.3f'%runtime,'s')
    return TProfile, cache



if __name__ == "__main__":
    para = parameter.main()
    results, cache = solve(para)
    

