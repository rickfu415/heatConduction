# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 23:07:48 2019

@author: RickFu
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd



def evolutionField(results):
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
    ax = fig.gca(projection='3d')
    ax.set_xlabel('x, m')
    ax.set_ylabel('Time, s')
    ax.set_zlabel('Temperature, K')
    Z = results.T.values
    ax.plot_surface(X, Y, Z, 
                    cmap=cm.seismic,
                    linewidth=0, 
                    antialiased=False)
    plt.show()



def thermalCouplePlot(results, positions):
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



def temperatureDistribution(results, times):
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
    
    length = parameter['length']
    numberOfNode = parameter['numberOfNode']
    numOfTimeStep = parameter['numberOfTimeStep']
    deltaTime = parameter['deltaTime']
    time = deltaTime * numOfTimeStep
    grids = np.linspace(0, length, numberOfNode).round(5)
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
    
    
    
    
    
    
    
    
    
    
    
    
    