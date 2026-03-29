# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:01:18 2019

@author: RickFu
"""
import numpy as np


def fixedValue(value, U2):
    """  Dirichlet boundary condition

    Assume that value of variable at BC is fixed.
    Please see any numerical analysis text book for details.
    
    Return: float
    """
    
    Ug = 2 * value - U2
    return float(np.asarray(Ug).flat[0])



def fixedGradient(q, k, dx, U1):
    """  Neumann boundary condition
    
    Assume that the resulted gradient at BC is fixed.
    Please see any numerical analysis text book for details.
    
    Return: float
    """
    
    Ug = q / k * 2 * dx + U1
    return float(np.asarray(Ug).flat[0])



def secondOrder(U, dx, Ug1, Ug2):
    """ Calculate second order derivative
    
    Centered differencing approximation.
    D2U/Dx2 = (U[i-1] - 2U[i] + U[i+1])/dx**2
    
    For BC nodes, use the values on ghost nodes.
    
    Ug1: value on ghost node at x=0
    Ug2: value on ghost node at x=L
    
    Please see any numerical analysis text book for details.
    
    Return: numpy array
    """
    
    d2U = np.zeros((U.size, 1))
    for i in range(0, U.size):
        if i==0:
            d2U[i] = (Ug1 - 2*U[i] + U[i+1]) / dx**2
        elif i==(U.size - 1):
            d2U[i] = (U[i-1] - 2*U[i] + Ug2) / dx**2
        else:
            d2U[i] = (U[i+1] - 2*U[i] + U[i-1]) / dx**2
    return d2U


def variableCoefficientDiffusion(U, k, dx_arr, Ug1, Ug2):
    """ Compute d/dx[k(x) * dU/dx] with variable conductivity and non-uniform grid.

    Vectorized implementation — no Python loops.

    dx_arr: array of length N-1, dx_arr[i] = distance from node i to node i+1.
    Ghost spacing mirrors the adjacent interior spacing.

    Uses arithmetic-mean interface conductivities:
        k_{i+1/2} = (k[i] + k[i+1]) / 2

    Return: numpy array (N, 1)
    """
    Uf = U.ravel()
    N = Uf.size

    # Interface conductivities: k_{i+1/2} for i=0..N-2, plus ghost interfaces
    k_east = np.empty(N)
    k_east[:-1] = 0.5 * (k[:-1] + k[1:])
    k_east[-1] = k[-1]  # ghost

    k_west = np.empty(N)
    k_west[1:] = 0.5 * (k[:-1] + k[1:])
    k_west[0] = k[0]  # ghost

    # Interface spacings
    dx_east = np.empty(N)
    dx_east[:-1] = dx_arr
    dx_east[-1] = dx_arr[-1]  # ghost mirror

    dx_west = np.empty(N)
    dx_west[1:] = dx_arr
    dx_west[0] = dx_arr[0]  # ghost mirror

    # Neighbor temperatures (with ghost values at boundaries)
    U_east = np.empty(N)
    U_east[:-1] = Uf[1:]
    U_east[-1] = Ug2

    U_west = np.empty(N)
    U_west[1:] = Uf[:-1]
    U_west[0] = Ug1

    # Flux balance
    h = 0.5 * (dx_west + dx_east)
    flux_east = k_east * (U_east - Uf) / dx_east
    flux_west = k_west * (Uf - U_west) / dx_west
    D = ((flux_east - flux_west) / h).reshape(-1, 1)
    return D



