# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 23:07:48 2019

@author: RickFu
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import os



def evolutionField(results, outputDir=None):
    """ Generate 3D temperature fields

    For better understanding of the results

    Inputs:
        1. parameter, a pandas series
        2. results, a numpy array
    """

    X = results.index
    Y = results.columns
    X, Y = np.meshgrid(X, Y)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('x, m')
    ax.set_ylabel('Time, s')
    ax.set_zlabel('Temperature, K')
    Z = results.T.values
    ax.plot_surface(X, Y, Z,
                    cmap=cm.seismic,
                    linewidth=0,
                    antialiased=False)
    if outputDir:
        fig.savefig(os.path.join(outputDir, 'evolutionField.png'), dpi=150)
    plt.close(fig)



def thermalCouplePlot(results, positions, outputDir=None):
    """ Generate x-y plots as thermo-couple data

    Inputs:
        1. results, a pandas DataFrame
        2. Positions, a list of positions of the generated
           grids.

    """

    df = results.loc[positions,:]
    df = df.T
    df = df.add_prefix('x = ')
    df = df.add_suffix(' m')
    ax = df.plot(grid=True)
    ax.set_xlabel("Time, s")
    ax.set_ylabel("Temperature, K")
    if outputDir:
        ax.get_figure().savefig(os.path.join(outputDir, 'thermalCouple.png'), dpi=150)
    plt.close(ax.get_figure())



def temperatureDistribution(results, times, outputDir=None):
    """ Generate temperature distribution at different times

    Inputs:
        1. results, a pandas DataFrame
        2. times, a list of timings on the calculated
           time steps
    """

    df = results.loc[:,times]
    df = df.add_prefix('t = ')
    df = df.add_suffix(' s')
    ax = df.plot(grid=True)
    ax.set_xlabel("x, m")
    ax.set_ylabel("Temperature, K")
    if outputDir:
        ax.get_figure().savefig(os.path.join(outputDir, 'temperatureDistribution.png'), dpi=150)
    plt.close(ax.get_figure())



def probePositions(parameter, probes_per_layer=5):
    """ Generate thermocouple probe positions from layer configuration.

    Places probes evenly within each layer, snapped to actual grid nodes.
    """
    # Build actual grid positions
    dx_arr = parameter['dx_array']
    N = parameter['numberOfNode']
    x = np.zeros(N)
    for i in range(N - 1):
        x[i + 1] = x[i] + dx_arr[i]

    def snap(target):
        """Snap target position to nearest grid node."""
        idx = np.argmin(np.abs(x - target))
        return round(x[idx], 5)

    if parameter.get('material function') == 'layered':
        t_layers = parameter['layerThicknesses']
        boundaries = np.cumsum(t_layers)
        starts = np.concatenate([[0], boundaries[:-1]])
        positions = set()
        for s, e in zip(starts, boundaries):
            for p in np.linspace(s, e, probes_per_layer):
                positions.add(snap(p))
        return sorted(positions)
    else:
        length = x[-1]
        positions = set()
        for p in np.linspace(0, length, probes_per_layer * 2):
            positions.add(snap(p))
        return sorted(positions)


def preprocess(parameter, results):
    """ Pre-Process results
    
    To convert numpy array into pandas DataFrame for easier
    data processing.
    
    Input:
        1. Generated parameter serie
        2. results as a numpy array
    
    Return:
        A pandas DataFrame with index as times and 
        columns as grid positions
    """
    
    numberOfNode = parameter['numberOfNode']
    numOfTimeStep = parameter['numberOfTimeStep']
    deltaTime = parameter['deltaTime']
    time = deltaTime * numOfTimeStep
    # Build x-positions from dx_array (supports non-uniform layered grid)
    dx_arr = parameter['dx_array']
    grids = np.zeros(numberOfNode)
    for i in range(numberOfNode - 1):
        grids[i + 1] = grids[i] + dx_arr[i]
    grids = grids.round(5)
    times = np.linspace(0, time, numOfTimeStep+1).round(5)
    df = pd.DataFrame(results, 
                      index = grids, 
                      columns = times)
    return df



if __name__ == "__main__":
    global para, results
    test = preprocess(para, results)
    evolutionField(test)
    positions = [0, 0.002, 0.004, 0.006, 0.008, 0.01]
    thermalCouplePlot(test, positions)
    times = [0, 2, 4, 6, 8, 10]
    temperatureDistribution(test, times)
    
    
    
    
    
    
    
    
    
    
    
    
    