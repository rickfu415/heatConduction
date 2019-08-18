# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 10:49:25 2019

@author: RickFu
"""
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import postprocessing as pp
import numpy as np
import parameter


def calRhs(para, x, t):
    """ Formula for T(x, t)
    
    Ref. Adam Joseph Amar, Modeling of One-Dimensional Ablation with Porous 
    Flow using Finite Control Volume Procedure.
    
    For details about the analytic solution, please read 
    Section 8.1.2, Eqs 8.9
    """
    
    length = para['length']
    conductivity = para['conductivity']
    density = para['density']
    heatCapacity = para['heatCapacity']
    alpha = conductivity / density / heatCapacity
    M_PI = 3.14159265358979323846
    TOL = 1E-18
    atol2 = alpha * t / length / length
    xol = x / length
    rhs = atol2 + 1. / 3 - xol + 0.5 * xol**2
    n = 1
    nPi = TOL + 1
    summ = 1
    while (summ > TOL):
        nPi = n * M_PI;
        summ = 2 / nPi**2 * np.exp(-nPi**2 * atol2) * np.cos(nPi * xol)
        rhs -= summ
        n += 1
    return rhs



def calT(para, x, t):
    """ Formula for T(x, t)
    
    Ref. Adam Joseph Amar, Modeling of One-Dimensional Ablation with Porous 
    Flow using Finite Control Volume Procedure.
    
    For details about the analytic solution, please read 
    Section 8.1.2, Eqs 8.9
    """
    
    T0 = para['IC value']
    qDotX0 = para['x=0 value']
    k = para['conductivity']
    length = para['length']
    T = T0 + qDotX0 * length * calRhs(para, x, t) / k
    return T



def solve(para):
    """ Obtain analytic solution
    
    Notes:
        1. Only work for heat flux boundary conditions
        2. Return the temperature profile
        3. Input para is the same parameters as used in 
           heatConduction solver
    
    Return: Numpy array containing temperature profile
    """
    
    numOfTimeStep = para['numberOfTimeStep']
    numberOfNode = para['numberOfNode']
    Tic = para['IC value']
    deltaTime = para['deltaTime']
    length = para['length']
    dx = length / numberOfNode
    grid = np.arange(0, length+dx, dx)
    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    T = np.zeros([numberOfNode, 1])
    TProfile[:,0] = Tic
    for timeStep in range(1, numOfTimeStep+1):
        print('[Analytic] Timestep',timeStep)
        time = timeStep * deltaTime
        for x in range(0, numberOfNode):
            T[x] = calT(para, grid[x], time)
        TProfile[:,timeStep] = T.reshape(1,-1)
    return TProfile



def compareError(times, T_numerical, T_analytic):
    """ Calculate L2 error
    
    Compare L2 error for all temperature distribution
    at given times
    
    Return a DataFrame
    """ 
    
    T_numerical = T_numerical[times]
    T_analytic = T_analytic[times]
    df = ((T_numerical - T_analytic)**2 / T_analytic**2)**0.5
    df = df.add_prefix('t = ')
    df = df.add_suffix(' s')
    ax = df.plot(grid=True)
    ax.set_xlabel("x, m")
    ax.set_ylabel("Relative error, %")
    return df



def plotDistribution(times, T_numerical, T_analytic):
    """ Plot and compare distribution
    
    Given times, numerical and analytic solutions, plot 2D 
    distributions
    """
    
    # Select results for given times
    T_numerical = T_numerical[times]
    T_analytic = T_analytic[times]
    
    # Obtain positions for text on figure
    x = T_numerical.index[5]
    y = T_numerical.loc[x,:].values
    
    # Adding column names
    T_numerical = T_numerical.add_prefix('t = ')
    T_numerical = T_numerical.add_suffix(' s')
    T_analytic = T_analytic.add_prefix('t = ')
    T_analytic = T_analytic.add_suffix(' s')
    
    # Plot numerical solutions
    ax = T_numerical.plot(color = 'r',
                          grid=True,
                          legend=False)
    
    # Adding text informations for each line
    cols = T_numerical.columns
    for i, col in enumerate(cols):
        t = ax.text(x,y[i],col)
        t.set_bbox(dict(alpha=0.5, 
                        edgecolor='red'))
    # Set axis names
    ax.set_xlabel("x, m")
    ax.set_ylabel("Temperature, K")
    
    # Plot analytic solutinos
    T_analytic.plot(ax=ax, marker='.', lw=0.1, 
            color = 'b',
            grid=True,
            legend=False)
    
    # Manually add legend
    legend_elements = [Line2D([0],[0], color='r', label='Nnumerical'),
                       Line2D([0],[0], marker='.', color='b', 
                              label='Analytic')]
    ax.legend(handles=legend_elements)




if __name__ == "__main__":
    global T
    para = parameter.main()
    analytics = solve(para)
    analytics = pp.preprocess(para, analytics)
    times = [2, 4, 6, 8, 10]
    plotDistribution(times, T, analytics)
    l2error = compareError(times, T, analytics)
    





