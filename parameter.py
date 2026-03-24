# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2019

@author: RickFu
"""
import postprocessing as pp
import heatConduction as hc
import pandas as pd
import numpy as np
import json
import os


MATERIAL_DB_DIR = os.path.join(os.path.dirname(__file__), 'material_db', 'simplified')


def load_material(name):
    """ Load material properties from JSON database.

    Input: material id string (e.g. 'pica', 'steel_304', 'li900')
    Return: dict with density, conductivity, heat_capacity
    """
    path = os.path.join(MATERIAL_DB_DIR, name + '.json')
    with open(path, 'r') as f:
        return json.load(f)


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

    # Material layers: list of material names from the database.
    # One name = constant material, multiple = layered TPS.
    df.at['materials'] = ['pica', 'steel_304']
    df.at['layerThicknesses'] = np.array([0.05, 0.05])  # sum = length

    # Derived from materials list (populated by load_materials)
    df = load_materials(df)
    
    # Grid: nodes per layer (boundaries are shared between layers)
    df.at['length'] = 0.1
    df.at['nodesPerLayer'] = 51  # nodes in each layer (including shared boundaries)
    
    # Solution
    df.at['numberOfTimeStep'] = 400#400
    df.at['deltaTime'] = 1
    df.at['maxIteration'] = 100
    df.at['convergence'] = 1E-9
    df.at['relaxation'] = 0.9 # value in [0-1] Very sensitive!!!
    
    # Environment
    df.at['ambientTemperature'] = 298.
    df.at['stefanBoltzmann'] = 5.670374419e-8  # W/(m^2 K^4)
    df.at['emissivity'] = 0.9

    # Initial conditions
    df.at['IC value'] = 298.

    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux'#'heatFlux' or 'fixedTemperature'
    df.at['x=0 value'] = 1e6
    df.at['x=L type'] = 'heatFlux'#'heatFlux' or 'fixedTemperature'
    df.at['x=L value'] = 0.
    df.at['re-radiate'] = True
    
    # Differential target
    df.at['back_wall_temperature_target'] = 450

    # Optimizer
    df.at['learning_rate'] = 1e-9  # thickness gradients are O(1e5), thicknesses O(1e-3)
    return df


def load_materials(df):
    """ Read material properties from JSON database based on materials list.

    Sets:
        material function: 'constant' or 'layered'
        numberOfLayers, layerConductivities (from DB)
        density, conductivity, heatCapacity (scalar for constant,
            or representative values for layered)
        layerDensities, layerHeatCapacities (arrays for layered)
    """
    materials = df.at['materials']
    n = len(materials)
    props = [load_material(m) for m in materials]

    if n == 1:
        df.at['material function'] = 'constant'
        df.at['density'] = props[0]['density']
        df.at['conductivity'] = props[0]['conductivity']
        df.at['heatCapacity'] = props[0]['heat_capacity']
    else:
        df.at['material function'] = 'layered'
        df.at['numberOfLayers'] = n
        df.at['layerConductivities'] = np.array([p['conductivity'] for p in props])
        df.at['layerDensities'] = np.array([p['density'] for p in props])
        df.at['layerHeatCapacities'] = np.array([p['heat_capacity'] for p in props])
        # Use first layer values as representative for solver parameters
        # that expect scalars (density, heatCapacity are overridden per-node
        # in normalize_conductivity)
        df.at['density'] = props[0]['density']
        df.at['conductivity'] = props[0]['conductivity']
        df.at['heatCapacity'] = props[0]['heat_capacity']
    return df


def normalize_conductivity(para):
    """ Build per-node property arrays from layer parameters.

    Layered grid: each layer gets its own uniform sub-grid with
    nodesPerLayer nodes. Adjacent layers share their boundary node.
    Total nodes = nodesPerLayer * n_layers - (n_layers - 1).

    When the optimizer changes layer thicknesses, dx within each
    layer adjusts — no nodes cross material boundaries.

    For constant material: broadcasts scalar to array.
    """
    if para['material function'] == 'layered':
        k_layers = para['layerConductivities']
        rho_layers = para['layerDensities']
        cp_layers = para['layerHeatCapacities']
        t_layers = para['layerThicknesses']
        n_layers = len(k_layers)
        npl = para['nodesPerLayer']

        # Total nodes: shared boundary nodes between layers
        N = npl * n_layers - (n_layers - 1)
        para['numberOfNode'] = N

        k_array = np.zeros(N)
        rho_array = np.zeros(N)
        cp_array = np.zeros(N)
        dx_array = np.zeros(N - 1)  # dx between consecutive nodes

        idx = 0
        for l in range(n_layers):
            n_nodes = npl
            dx_l = t_layers[l] / (n_nodes - 1)
            # Fill properties for this layer's nodes
            for i in range(n_nodes):
                node = idx + i
                k_array[node] = k_layers[l]
                rho_array[node] = rho_layers[l]
                cp_array[node] = cp_layers[l]
            # Fill dx for this layer's intervals
            for i in range(n_nodes - 1):
                dx_array[idx + i] = dx_l
            # Next layer starts at the shared boundary node
            idx += n_nodes - 1

        para['conductivity'] = k_array
        para['density'] = rho_array
        para['heatCapacity'] = cp_array
        para['dx_array'] = dx_array
    else:
        N = para['numberOfNode']
        k = para['conductivity']
        if np.isscalar(k):
            para['conductivity'] = np.full(N, k)
        rho = para['density']
        if np.isscalar(rho):
            para['density'] = np.full(N, rho)
        cp = para['heatCapacity']
        if np.isscalar(cp):
            para['heatCapacity'] = np.full(N, cp)
        dx = para['length'] / (N - 1)
        para['dx_array'] = np.full(N - 1, dx)
    return para



if __name__ == "__main__":
    parameter = main()
    outputDir = parameter['output']
    os.makedirs(outputDir, exist_ok=True)
    results, cache = hc.solve(parameter)
    T = pp.preprocess(parameter, results)
    T.to_csv(os.path.join(outputDir, 'solutionHistory.csv'))
    cache['Log'].to_csv(os.path.join(outputDir, 'solverLog.csv'))
    pp.evolutionField(T, outputDir)
    positions = pp.probePositions(parameter)
    pp.thermalCouplePlot(T, positions, outputDir)
    totalTime = parameter['numberOfTimeStep'] * parameter['deltaTime']
    times = np.linspace(0, totalTime, 6).round(5).tolist()
    pp.temperatureDistribution(T, times, outputDir)
    
    
    
    
    