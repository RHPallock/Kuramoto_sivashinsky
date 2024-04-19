# Rumayel Hassan Pallock
# Kuramoto sivanski Function and linear operator
# Date 11/22/2023
# equation----- dt(u) = -u*dx(u) - dx(dx(u)) - dx(dx(dx(dx(u))))


from tracemalloc import stop
import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import h5py # HDF5 file manipulation
from mpi4py import MPI
import logging
logger = logging.getLogger(__name__)
import shutil
shutil.rmtree('full_solution',ignore_errors = True)

#Parameters
L = 22

# Define co-ordinate, distributor and basis
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord,dtype = np.float64)
xbasis = d3.RealFourier(xcoord,size = 32,bounds = (0,L),dealias = 1)


# RHS operator of KS equation
def RHS(u_b):

    ux = d3.Differentiate(u_b,xcoord)
    uxx = d3.Differentiate(ux,xcoord)
    uxxx = d3.Differentiate(uxx,xcoord)
    uxxxx = d3.Differentiate(uxxx,xcoord)

    F_u = (-u_b*ux - uxx - uxxxx)
    return F_u
def RHS_ev(u_b):

    ux = d3.Differentiate(u_b,xcoord)
    uxx = d3.Differentiate(ux,xcoord)
    uxxx = d3.Differentiate(uxx,xcoord)
    uxxxx = d3.Differentiate(uxxx,xcoord)

    F_u = (-u_b*ux - uxx - uxxxx).evaluate()
    return F_u

def Adjoint_Operator(u_b,u_pp):

    u_pp_x = d3.Differentiate(u_pp,xcoord)
    u_pp_xx = d3.Differentiate(u_pp_x,xcoord)
    u_pp_xxx = d3.Differentiate(u_pp_xx,xcoord)
    u_pp_xxxx = d3.Differentiate(u_pp_xxx,xcoord)

    Adj_L = (u_b*u_pp_x - u_pp_xx - u_pp_xxxx)
    return Adj_L



# define field variable u
u = dist.Field (name='u', bases = xbasis)
# define the derivatives with substituions
dx = lambda A: d3.Differentiate(A,xcoord)
F_u = lambda A: RHS(A)
Adj = lambda A,B: Adjoint_Operator(A,B)

# Add equations
# main equation , with linear terms on lhs and nonlinear terms on rhs
problem = d3.IVP([u],namespace = locals())
problem.add_equation("dt(u)+dx(dx(dx(dx(u))))+2*dx(dx(dx(dx(dx(dx(u)))))) + dx(dx(dx(dx(dx(dx(dx(dx(u)))))))) = u*dx(u*dx(u)+dx(dx(u))+dx(dx(dx(dx(u)))))-dx(dx(u*dx(u))) - dx(dx(dx(dx(u*dx(u)))))")

# Initial Guess
x = dist.local_grid(xbasis)
u['g'] = 2*np.sin(2*np.pi*x/L)

# Build solver
solver = problem.build_solver(d3.RK443)

# Setting stop criteria for fictitous time
stop_time = 4000
solver.stop_sim_time =stop_time

# analysis
# Output solution fields
full_solution = solver.evaluator.add_file_handler('full_solution', iter =25 , max_writes = 500000)
full_solution.add_task(u, layout = 'g', name='u')




# main loop
timestep = .1
total_iteration = int(stop_time/timestep)
L2_residue = np.zeros((total_iteration,1))
iteration = 0
while solver.proceed:
    solver.step(timestep)
    logger.info('L2 residue = %.10f' %np.linalg.norm([RHS_ev(u)['g']*RHS_ev(u)['g']]))
    L2_residue[iteration,0] = np.linalg.norm([RHS_ev(u)['g']*RHS_ev(u)['g']])
    iteration = iteration+1
    if solver.iteration % 100 == 0:
        logger.info('Iteration = %i, time = %e, dt = %e' %(solver.iteration, solver.sim_time, timestep))
    if L2_residue[iteration-1,0]<=1e-10:
        break


file = h5py.File('L2_Residue_E3.h5','w')
file.create_dataset('L2_Residue',data = L2_residue)
file.close()

u.change_scales(1)
print(u['g'])
plt.plot(x,u['g'])
plt.show()

fic_time = np.linspace(0,total_iteration-1,total_iteration)
plt.plot(fic_time,L2_residue[0:total_iteration,0])
plt.xlabel("Iteration")
plt.ylabel("L2 residue")
plt.show()
