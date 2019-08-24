# 1D Heat Conduction Solver (08-21-2019)
A transient 1D heat conduction solver using Finite Difference Method and implicit backward Euler time scheme.  

## Features:  
    1. Fully modularized, easy to customize for your own problem.  
    1. Only use the common packages, Numpy, Pandas and Matplotlib.  
    2. Centered Differecing in space (second order accuracy), implicit backward Euler time scheme (First order accuracy).  
    3. Using Newton's method to solve discretized equation system at each time step.  
    4. Two types of boundary conditions, fixed temperature and fixed heat flux.  
    5. Current version only support constant material properties, will be upgraded.  
    6. Incoporate with postprocessing module and analytic solution comparison.   

## How to run:  
    1. In any Python IDE, open parameter.py, execute.  
    2. To compare with analytic solution, open analyticSolution.py, execute.  
    
## Reference:
1. [Numerical analysis](http://web.engr.uky.edu/~acfd/egr537-lctrs.pdf).
2. [Modeling of One-Dimensional Ablation with Porous Flow using Finite Control Volume Procedure](http://www.lib.ncsu.edu/resolver/1840.16/2847).

## Citation:
If you are using the code for your research or study, please consider cite me :) I am a miserable PhD ......
1. [Thermomechanical Coupling for Charring Ablators](https://doi.org/10.2514/1.T5194).
2. [Thermal Expansion for Charring Ablative Materials](https://doi.org/10.2514/1.T5718).
