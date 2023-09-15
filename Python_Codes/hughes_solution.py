import sys
import math
from math import *
from dolfin import *
import ufl
from ufl import ge, eq
from numpy import array, heaviside

# Set some fenics parameters
# parameters["form_compiler"]["quadrature_degree"] = 150
# dolfin.parameters['ghost_mode'] = 'shared_facet'

from setup_example import *
##############################################################
#     Define options, model and discretization parameters    #
##############################################################

# Define options
with_plot = True               # Enable plotting
verbosity = 2                  # Verbosity: 1=less, 2=more output

# Define directory for plots and log files
output_basedir = 'output'

# Choose example (defined in setup_example.py)
example = Example3D1()

# Setup model parameters
delta_1 = Constant(2.e-1)      # Laplacian in Eikonal eq.
delta_2 = Constant(1.e-1)      # Right-hand side of Eikonal eq.

epsilon = Constant(1.e-5)      # Laplacian in rho eq.
penalty = Constant(1.e0)       # jump penalty term in DG formulation (must be one for zero-order)

outflow_vel = 10.0             # Outflow velocity (gamma in the boundary condition)

# Note: More parameters are defined in the example class (see setup_example.py)

##############################################################
#       Main program starts here (don't cross this line)     #
##############################################################
# Define parameters for time discretization
T           = example.T                          # Final time T
N           = example.N                        # Number of time steps
tau         = T/N                      # Time discretization parameter

# Get data from example
mesh        = example.mesh


rho_0       = example.rho_0
ag_pos_0    = example.ag_pos_0


BoundingBox = mesh.bounding_box_tree()

# Indicator functions for the boundaries
boundary = lambda x,on_boundary : example.exits.inside(x, on_boundary)

nr_agents = np.size(ag_pos_0, 0)

# 3D
cell = ufl.Cell("triangle", 3) # triangle cell, embedded in 3D space
coordinate_element_rho = ufl.FiniteElement("DG", cell, 0)
coordinate_element_phi = ufl.FiniteElement("CG", cell, 1)
R = FunctionSpace(mesh, coordinate_element_rho)           # Density
P = FunctionSpace(mesh, coordinate_element_phi)           # Potential


# Dirichlet boundary condition for phi
bc_phi = DirichletBC(P, Constant(0), boundary)

# Set some trial and test functions
phi_ = TrialFunction(P)
z = TestFunction(P)
rho_ = TrialFunction(R)
w = TestFunction(R)

# Mesh entities
n = FacetNormal(mesh)
h = CellDiameter(mesh)
x = SpatialCoordinate(mesh)

# Define new surface measure for the boundary conditions
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

boundary_markers.set_all(0)

exits = example.exits
exits.mark(boundary_markers, 1)

dss = Measure('ds', domain=mesh, subdomain_data=boundary_markers) # Exits

sip = Constant(1.0)

# Get UFL expression of the numerical flux
def Flux(rho_, beta, n):
    
    # Lax-Friedrichs flux
    flux = inner(avg(rho_*beta), n('+')) - 0.5*jump(rho_)
    return -flux
  
# Potential function
def f(rho):
    return 1.-rho
def Df(rho):
    return -1.0

# Velocity constraint
def v_trunc(Dphi):
    return conditional(lt(inner(Dphi,Dphi),1), Dphi, Dphi/sqrt(inner(Dphi,Dphi)))


def SolveForward():

    rhos = np.empty((N+1,len(R.dofmap().dofs())))
    phis = np.empty((N+1,len(P.dofmap().dofs())))
    
    rho = interpolate(rho_0, R)
    rhos[0,:] = rho.vector().get_local()

    rho_old = Function(R)
    rho_old.assign(rho)

    # Compute initial potential
    phi = Function(P)    
    phi_ = TrialFunction(P)    
        
    # Define final phi equation
    eq_phi = delta_1 * inner(grad(phi), grad(z))*dx \
             + inner(grad(phi), grad(phi))*z*dx \
             - (1.0/(f(rho)**2 + delta_2))*z*dx

    # Solve Eikonal equation    
    J  = derivative(eq_phi, phi, phi_)   # Gateaux derivative in dir. of phi_
    problem = NonlinearVariationalProblem(eq_phi, phi, bc_phi,J)
    solver  = NonlinearVariationalSolver(problem)

    prm = solver.parameters
    prm['nonlinear_solver'] = 'snes'
    prm['snes_solver']['line_search'] = 'bt'
    prm['snes_solver']['linear_solver'] = 'mumps'
    prm['snes_solver']['preconditioner'] = 'hypre_amg'
    prm['snes_solver']['krylov_solver']['nonzero_initial_guess'] = False
    prm['snes_solver']['report'] = True

    solver.solve()
    
    phis[0,:] = phi.vector().get_local()
    ############################
    # Solve the state equation #
    ############################
    if verbosity > 1:
        print("  solving forward equation")

    for k in range(N):
        
        if verbosity > 1:
            print(" time step ", k+1, " of ", N, end="\r")
        #######################################
        # Solve convection diffusion equation #
        #######################################
        # Velocity term
        vv = Expression(("1","-1","0"),  domain=mesh,degree=2)

	#Another velocity for controlling the direction
        v1 = Expression(("-5*exp(-pow(x[2]-1,2))+5*exp(-pow(x[2]+1,2))",'0',"-10*exp(-pow(x[0]+5,2))"),  domain=mesh,degree=2)
        v2 = Expression(("-5*exp(-pow(x[2]-1,2))",'0',"-10*exp(-pow(x[0]+5,2))"),  domain=mesh,degree=2)
        vv1 = conditional(lt(x[0]-3.75,0),v1,v2)

        # Defining Convection terms regarding vv and vv1
        conv_term = inner(vv,avg(w)*jump(rho_,-n))*dS 
        conv_term1 = inner(vv1,avg(w)*jump(rho_,-n))*dS 

        
        # Bilinear form for DIFFUSION
             # dx ... domain integral
             # dS ... interior facets integral
             # ds ... exterior facets integral
             # dss ... custom measure, dss(1)
        a_sip = epsilon * inner(grad(rho_), grad(w))*dx\
                - sip*epsilon * inner(avg(grad(rho_)),jump(w,n))*dS \
                - sip*epsilon * inner(avg(grad(w)),jump(rho_,n))*dS \
                + sip*(penalty/avg(h)) * epsilon * jump(rho_)*jump(w)*dS \
                + sip*outflow_vel*rho_*w*dss(1) \
                + conv_term1

        
        # Transport direction
        beta = f(rho_old)*v_trunc(grad(phi))
       
        # Bilinear form for TRANSPORT
        flux = Flux(rho_old, beta, n)
   
        a_upw = flux*jump(w)*dS(metadata = {"quadrature_degree": 5})
        
        # Bilinear form for TIME DISCRETIZATION
        a_time = rho_*w*dx
        
        
        a_full = a_time + tau*a_sip 
        L_rho = rho_old*w*dx -tau*a_upw 
        
        # Solve system and store solution
        solve(a_full== L_rho, rho)
        rhos[k+1,:] = rho.vector().get_local()
        rho_old.assign(rho)
        ##############################
        # Solve the Eikonal equation #
        ##############################
        if k <= N-1:
           solve(eq_phi == 0, phi, bc_phi)
           phis[k+1,:] = phi.vector().get_local()
           phi.assign(phi)
        sys.stdout.flush()

    if verbosity > 1:
        print('\n')
        
    return [rhos, phis]

