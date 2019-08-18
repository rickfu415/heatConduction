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



