# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2019

@author: RickFu
"""
import postprocessing as pp
import heatConduction as hc
import pandas as pd
import os


def main():
    """ Generate parameter
    
    1. Generate system-level parameters
    2. Generate material properties, grid, time, bcs
    
    Return: a Pandas series
    """
    
    column = 'values'
    df = pd.Series(name = column)
    df = df.astype('object')
    
    # System-level 
    df.at['problem'] = 'HeatConduction'
    df.at['SpatialDiscretize'] = 'CenteredDifferencing'
    df.at['TimeDiscretize'] = 'BackwardEular'
    df.at['ODEsolver'] = 'NewtonIteration'
    df.at['linearSolver'] = 'numpy linalg'
    df.at['CPU'] = 1
    df.at['output'] = 'results'
    
    # Material
    df.at['material'] = 'steel'
    df.at['material function'] = 'constant'
    df.at['density'] = 7850
    df.at['conductivity'] = 60.5
    df.at['heatCapacity'] = 434
    
    # Grid
    df.at['length'] = 0.01
    df.at['numberOfNode'] = 101
    
    # Solution
    df.at['numberOfTimeStep'] = 100#400
    df.at['deltaTime'] = 0.1
    df.at['maxIteration'] = 100
    df.at['convergence'] = 1E-9
    df.at['relaxation'] = 0.9 # value in [0-1] Very sensitive!!!
    
    # Initial conditions
    df.at['IC value'] = 298.
    
    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux'#'heatFlux' or 'fixedTemperature'
    df.at['x=0 value'] = 750000
    df.at['x=L type'] = 'heatFlux'#'heatFlux' or 'fixedTemperature'
    df.at['x=L value'] = 0.

    # Differential target
    df.at['back_wall_temperature_target'] = 450
    return df



if __name__ == "__main__":
    parameter = main()
    outputDir = parameter['output']
    os.makedirs(outputDir, exist_ok=True)
    results, cache = hc.solve(parameter)
    T = pp.preprocess(parameter, results)
    T.to_csv(os.path.join(outputDir, 'solutionHistory.csv'))
    cache['Log'].to_csv(os.path.join(outputDir, 'solverLog.csv'))
    pp.evolutionField(T, outputDir)
    positions = [0, 0.002, 0.004, 0.006, 0.008, 0.01]
    pp.thermalCouplePlot(T, positions, outputDir)
    times = [0, 2, 4, 6, 8, 10]
    pp.temperatureDistribution(T, times, outputDir)
    
    
    
    
    