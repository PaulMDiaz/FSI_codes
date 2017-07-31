#  Serves as the Driver for FSI analysis of Cavity
#  Driven Flow with a bottom non-rigid wall
#
#  Created: 26 February 2014
#  Modified: 26 July 2016

# Import user-defined solvers for fluid and structure domains
from StructureSolver import Structure_Solver
from FluidSolver import Fluid_Solver
from DrivenCavity import DrivenCavity

# Import IO object
from iotools import *

# Import dolfin and numpy
from fenics import *
import numpy as np
import pylab as plt
import math

parameters["allow_extrapolation"] = True
set_log_level(ERROR)

############## INPUT DATA PARAMETERS ###################

# Physical parameters
nu_f = 0.01	# Fluid viscosity
nu_s = 0.2	# Structure Poisson coefficient
E_s = 1e3	# Structure Young modulus
rho_f = 1	# Fluid density (incorporated in the fluid corrected pressure as p_corr = p/rho)

# Numerical parameters
dt = 0.025	# Time step
T = 1		#  Set final time for iteration
N = 32		# Number of discretizations (square mesh)

# Geometric parameters
h = 0.25	# Nondimensional structure thickness
# Check if N is a multiple of 1/h -> Error check to be included

#  Create Object with Driven Cavity expressions, mesh and subdomains
DC = DrivenCavity(nu_f, nu_s, E_s, rho_f, dt, T, N, h)
IO = IO()

############# Fluid Solver Setup section #################
#  Set the Fluid Element Type
FluidElementType = "Lagrange"

#  Set the VelocitySpaceDegree
VelocityElementDegree = 2

#  Set the PressureSpaceDegree
PressureElementDegree = 1

# Set the solver used for the fluid problem
# The only option available for the fluid solver is Chorin decomposition
FluidSolverMethod = "Chorin"

# Set the body forces on the fluid
FluidBodyForce = Constant((0.0,0.0))

#  Set the equations for the variational solver and boundary conditions and create fluid solver object
F = Fluid_Solver(DC.mesh_f, FluidElementType,
		VelocityElementDegree, PressureElementDegree, FluidSolverMethod, FluidBodyForce)

########### Structure Solver Setup section #################

#  Set the Structure Element Type
StructureElementType = "CG"

#  Set the Structure Element Degree
StructureElementDegree = 1

# Set the solver used for the structure problem
# The options are "Linear" or "NeoHookean" for the structure solver
StructureSolverMethod = "NeoHookean"
# Body forces on the structure
StructureBodyForce = Constant((0.0, 0.0))

# Set the equations for the variational solver and boundary conditions and creat structure solver object
S = Structure_Solver(DC.mesh_s, StructureElementType,
	StructureElementDegree, StructureSolverMethod, StructureBodyForce)

# Define boundary conditions on subdomains (fluid and structure) but traction from fluid on the structure which
# will be included after computing the fluid.
DC.Define_Boundary_Conditions(S, F)

################ Iteration section ##############################

# Sequentialy staggered iteration scheme
while DC.t < DC.T + DOLFIN_EPS:
	print ' STARTING TIME LOOP ITERATION ', DC.iter_t
	DC.t += DC.dt
	for ii in range(3):
		print ''
		print ''
		print 'Time loop iteration number = ', DC.iter_t
		print 'Loop iteration time = ', DC.t

		# Compute structural displacement and velocity
		S.Structure_Problem_Solver(DC, F)
		# Compute velocity mesh
		IO.Move_Mesh(S, F)
		# Solve fluid problem for velocity and pressure
		F.Fluid_Problem_Solver(DC, S)

	S.d0.assign(S.d)
	F.u0.assign(F.u1)
	for x in F.mesh.coordinates(): x[:] += DC.dt*F.u_mesh(x)[:]
	DC.Save_Results(S, F)
