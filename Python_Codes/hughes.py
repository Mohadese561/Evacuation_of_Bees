import sys, os
import time as timer
import math
import re
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from hughes_solution import *

# Advanced: uncomment to start with a previously computed control stored in file
# control_file = "output/control.csv"

set_log_level(LogLevel.WARNING)

np.set_printoptions(threshold=sys.maxsize)                        

# exit(0)
# File output
res_file = open(output_basedir + '/result.dat', "w+")

# Write evolution of objective to file
def WriteObjective(rhos):
        objective_file = open(output_basedir + '/objective.dat', 'w')
        
        for k in range(1,N+1):
                rho = Function(R)
                rho.vector().set_local(rhos[k,:])

                J_track = assemble(rho*dx(1))                

                objective_file.write(str(J_track) + "\n")

        objective_file.close()                

# Plot density
def plot_rho(rhos):
    
    if with_plot:
        file_rho = File(output_basedir + '/rho.pvd')
        file_mass = open(output_basedir + '/mass.csv', "w")
        
        for k in range(N):
            rho = Function(R)
            rho.vector().set_local(rhos[k,:])
            rho.rename("rho", "density")
            file_rho << rho,k

# Plot potential
def plot_phi(phis):
    
    if with_plot:
        file_phi = File(output_basedir + '/phi.pvd')

        for k in range(N):
            phi = Function(P)
            phi.vector().set_local(phis[k,:])

            phi.rename("phi", "potential")
            file_phi << phi,k

##############################
# MAIN program starts here
#############################     
# Solve forward equation
[rhos, phis] = SolveForward()

# Write current iterate to files
if verbosity > 1:
    print("  writing output files")

plot_rho(rhos)
plot_phi(phis)
#WriteObjective(rhos)
                
        
res_file.close()
