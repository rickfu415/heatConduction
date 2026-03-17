"""
Created on 03 17 2026

@author: RickFu
"""
import postprocessing as pp
import heatConduction as hc
import pandas as pd
import os
import parameter as parameter

def main(para, cache):
    """ 
    Main function


    """
    print("Start Adjoint Calculation")
    target_temperature = para['back_wall_temperature_target']






    return

if __name__ == "__main__":
    para = parameter.main()
    outputDir = para['output']
    print('Output directory: ' + outputDir)
    os.makedirs(outputDir, exist_ok=True)
    results, cache = hc.solve(para)
    main(para, cache)