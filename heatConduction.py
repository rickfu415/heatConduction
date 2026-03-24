# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:15:28 2019

@author: RickFu
"""
import numpy as np
import pandas as pd
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

    # Precompute interface conductivities and spacings
    k_half = np.zeros(numberOfNode + 1)
    dx_half = np.zeros(numberOfNode + 1)
    for i in range(numberOfNode - 1):
        k_half[i + 1] = 0.5 * (k[i] + k[i + 1])
        dx_half[i + 1] = dx_arr[i]
    k_half[0] = k[0];                dx_half[0] = dx_arr[0]
    k_half[numberOfNode] = k[-1];    dx_half[numberOfNode] = dx_arr[-1]

    # Compute effective heat flux at boundaries (including re-radiation)
    dx0 = dx_arr[0]
    dxN = dx_arr[-1]
    qX0 = valueX0
    qXL = valueXL
    if reradiate:
        if typeX0 == 'heatFlux':
            qX0 = valueX0 - eps_r * sigma * (T[0, 0]**4 - T_amb**4)
        if typeXL == 'heatFlux':
            qXL = valueXL - eps_r * sigma * (T[-1, 0]**4 - T_amb**4)

    # Boundary ghost values for temperature
    if typeX0 == 'heatFlux':
        Ug1 = utility.fixedGradient(qX0, k[0], dx0, T[1])
    elif typeX0 == 'fixedTemperature':
        Ug1 = utility.fixedValue(valueX0, T[1])

    if typeXL == 'heatFlux':
        Ug2 = utility.fixedGradient(qXL, k[-1], dxN, T[-2])
    elif typeXL == 'fixedTemperature':
        Ug2 = utility.fixedValue(valueXL, T[-2])

    # Assemble Jacobian node by node, store ce/cw for adjoint reuse
    ce_arr = np.zeros(numberOfNode)
    cw_arr = np.zeros(numberOfNode)
    for i in range(numberOfNode):
        k_east = k_half[i + 1]
        k_west = k_half[i]
        dx_e = dx_half[i + 1]
        dx_w = dx_half[i]
        h_i = 0.5 * (dx_w + dx_e)

        ce = dt / (rho[i] * hcp[i] * h_i * dx_e)
        cw = dt / (rho[i] * hcp[i] * h_i * dx_w)
        ce_arr[i] = ce
        cw_arr[i] = cw

        # Diagonal
        Jacobian[i][i] = 1 + ce * k_east + cw * k_west

        if i == 0:
            if typeX0 == 'heatFlux':
                Jacobian[0][1] = -(ce * k_east + cw * k_west)
            elif typeX0 == 'fixedTemperature':
                Jacobian[0][1] = -ce * k_east + cw * k_west
        elif i == numberOfNode - 1:
            if typeXL == 'heatFlux':
                Jacobian[-1][-2] = -(ce * k_east + cw * k_west)
            elif typeXL == 'fixedTemperature':
                Jacobian[-1][-2] = ce * k_east - cw * k_west
        else:
            Jacobian[i][i + 1] = -ce * k_east
            Jacobian[i][i - 1] = -cw * k_west

    # Re-radiation Jacobian correction at boundary nodes
    # dq_rad/dT = -4*eps*sigma*T^3, enters through the ghost flux
    # Extra diagonal contribution: dt/(rho*cp) * 2/(h_i) * eps*sigma*4*T^3
    if reradiate:
        if typeX0 == 'heatFlux':
            h_0 = 0.5 * (dx_half[0] + dx_half[1])
            rad_jac_0 = dt / (rho[0]*hcp[0]) * 2.0 / h_0 * eps_r * sigma * 4 * T[0, 0]**3
            Jacobian[0][0] += rad_jac_0
            Jacobian[0][1] -= rad_jac_0  # ghost depends on T[1] for heatFlux BC
        if typeXL == 'heatFlux':
            h_N = 0.5 * (dx_half[-2] + dx_half[-1])
            rad_jac_N = dt / (rho[-1]*hcp[-1]) * 2.0 / h_N * eps_r * sigma * 4 * T[-1, 0]**3
            Jacobian[-1][-1] += rad_jac_N
            Jacobian[-1][-2] -= rad_jac_N

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
    """ Solve Ax=B
    
    Process:
        1. Get A = Jacobian matrix (Jacobian)
        2. Get B = Right hand side equation (F)
        3. Calculate dT
        4. Update T
        5. Store in cache
        
    Return: a dictionary
    """
    relax = para['relaxation']
    A = cache['Jacobian']
    B = cache['F']
    dT = np.linalg.solve(A, B)
    T = cache['T']
    T = dT * relax + T
    cache['T']=T
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
    if verbose:
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
    if verbose:
        print(' [Step] [Pysical Time] [Iteration] [Residue]    [T_0]      [T_L]')
    for timeStep in range(1, numOfTimeStep+1):
        cache['ts'] = timeStep
        cache = newtonIteration(para, cache, verbose=verbose)
        cache = storeUpdateResult(cache)
    TProfile = cache['TProfile']
    runtime = time.time() - start
    if verbose:
        print('[Cost] CPU time spent','%.3f'%runtime,'s')
    return TProfile, cache



if __name__ == "__main__":
    para = parameter.main()
    results, cache = solve(para)
    

