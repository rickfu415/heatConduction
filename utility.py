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
    return Ug



def fixedGradient(q, k, dx, U1):
    """  Neumann boundary condition
    
    Assume that the resulted gradient at BC is fixed.
    Please see any numerical analysis text book for details.
    
    Return: float
    """
    
    Ug =  q / k * 2 * dx  + U1
    return Ug



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

    dx_arr: array of length N-1, dx_arr[i] = distance from node i to node i+1.
    Ghost spacing mirrors the adjacent interior spacing.

    Uses arithmetic-mean interface conductivities:
        k_{i+1/2} = (k[i] + k[i+1]) / 2

    At node i, the control volume has width h_i = (dx_west + dx_east) / 2,
    and the flux balance gives:
        D_i = [k_east*(U_{i+1}-U_i)/dx_east - k_west*(U_i-U_{i-1})/dx_west] / h_i

    Return: numpy array (N, 1)
    """
    N = U.size
    D = np.zeros((N, 1))
    for i in range(N):
        # East interface
        if i < N - 1:
            k_east = 0.5 * (k[i] + k[i + 1])
            dx_east = dx_arr[i]
        else:
            k_east = k[i]
            dx_east = dx_arr[-1]  # ghost mirror

        # West interface
        if i > 0:
            k_west = 0.5 * (k[i - 1] + k[i])
            dx_west = dx_arr[i - 1]
        else:
            k_west = k[i]
            dx_west = dx_arr[0]  # ghost mirror

        # Temperature neighbors
        if i == 0:
            U_west = Ug1
            U_east = U[1]
        elif i == N - 1:
            U_west = U[i - 1]
            U_east = Ug2
        else:
            U_west = U[i - 1]
            U_east = U[i + 1]

        h_i = 0.5 * (dx_west + dx_east)
        flux_east = k_east * (U_east - U[i]) / dx_east
        flux_west = k_west * (U[i] - U_west) / dx_west
        D[i] = (flux_east - flux_west) / h_i
    return D



