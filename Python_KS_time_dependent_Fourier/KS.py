# Rumayel Hassan Pallock
# Kuramoto sivanski equation solution using dedalus version 3.0
# Date 1/13/2023
# equation----- dt(u) = -u*dx(u) - dx(dx(u)) - dx(dx(dx(dx(u))))


import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import logging
import pathlib
import subprocess
import os
import h5py # HDF5 file manipulation
import scipy.io
from mpi4py import MPI
import time

import logging
logger = logging.getLogger(__name__)
#%matplotlib inline
#%config InlineBackened.figure_format = 'retina'

import shutil
shutil.rmtree('full_solution',ignore_errors = True)

#Parameters
L = 22


# Define co-ordinate, distributor and basis
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord,dtype = np.float64)
xbasis = d3.RealFourier(xcoord,size = 50,bounds = (-L/2,L/2),dealias = 3/2)


# define fields 
# define field variable u
u = dist.Field (name='u', bases = xbasis)
problem = d3.IVP([u],namespace = locals())

# define the derivatives with substituions
dx = lambda A: d3.Differentiate(A,xcoord)

# Add equations
# main equation , with linear terms on lhs and nonlinear terms on rhs
problem.add_equation("dt(u)+dx(dx(u))+dx(dx(dx(dx(u)))) = -u*dx(u)")


# Initial condition
x = dist.local_grid(xbasis)
u['g'] = np.arccosh(x)


# Build solver
solver = problem.build_solver(d3.RK222)

# Setting stop criteria
solver.stop_sim_time = 1000


# analysis
# Output solution fields
full_solution = solver.evaluator.add_file_handler('full_solution', iter =25 , max_writes = 500)
full_solution.add_task(u, layout = 'g', name='u')




# main loop
timestep = 0.01
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 100 == 0:
        logger.info('Iteration = %i, time = %e, dt = %e' %(solver.iteration, solver.sim_time, timestep))


