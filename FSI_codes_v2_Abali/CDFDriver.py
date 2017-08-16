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



#  Create Object with Driven Cavity expressions, mesh and subdomains
DC = DrivenCavity()
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
StructureElementType = "CG" # CG is equivilent to lagrange.

#  Set the Structure Element Degree
StructureElementDegree = 1 # Standard linear lagrange element.

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


### for printing the interaction

# Fluid  array of coordinates

# !!! have to check that these coordinates remain on boundary... don't move. Otherwise condition is wrong !!!

# S scalar, V vector, T Tensor
dofs_f_S = F.S_space.tabulate_dof_coordinates().reshape((F.S_space.dim(),-1))
dofs_f_V = F.V_space.tabulate_dof_coordinates().reshape((F.V_space.dim(),-1))
dofs_f_T = F.T_space.tabulate_dof_coordinates().reshape((F.T_space.dim(),-1))

# Structure array of coordinates
dofs_s_S = S.S_space.tabulate_dof_coordinates().reshape((S.S_space.dim(),-1))
dofs_s_V = S.V_space.tabulate_dof_coordinates().reshape((S.V_space.dim(),-1))
dofs_s_T = S.T_space.tabulate_dof_coordinates().reshape((S.T_space.dim(),-1))

#Extract dof indices for values on boundary.
# y = 0.5 if mesh is not deformed.
i_f_S = np.where((dofs_f_S[:,1] == 0.5))[0] #  & (x <= 0.5)
i_f_V = np.where((dofs_f_V[:,1] == 0.5))[0] #  & (x <= 0.5)
i_f_T = np.where((dofs_f_T[:,1] == 0.5))[0] #  & (x <= 0.5)

i_s_S = np.where((dofs_s_S[:,1] == 0.5))[0] #  & (x <= 0.5)
i_s_V = np.where((dofs_s_V[:,1] == 0.5))[0] #  & (x <= 0.5)
i_s_T = np.where((dofs_s_T[:,1] == 0.5))[0] #  & (x <= 0.5)

################ Iteration section ##############################
stop = 2*DC.dt
#print stop

# Sequentialy staggered iteration scheme
while DC.t < DC.T + DOLFIN_EPS:
	print ' STARTING TIME LOOP ITERATION ', DC.iter_t
	DC.t += DC.dt
	# Loop for convergence between fluid and structure
	for ii in range(3):
		print ''
		print ''
                # not sure what time loop iteration number is... how many to converge at given time step?
		print 'Time loop iteration number = ', DC.iter_t
		print 'Loop iteration time = ', DC.t

                # Solve fluid problem for velocity and pressure
		F.Fluid_Problem_Solver(DC, S)
		# Compute structural displacement and velocity
		S.Structure_Problem_Solver(DC, F)
		# Compute velocity mesh
		IO.Move_Mesh(S, F)

	S.d0.assign(S.d) # set current time to previous for displacement
	F.u0.assign(F.u1)# set current time to previous for velocity. (pressure?)

	if DC.t == stop:
		print "Checking interface solutions at time:", stop

		# Fluid velocity, pressure and stress
		u_FSI = F.u1.vector()[i_f_V]
        #u_FSI = F.v_.vector()[i_f_V]

		p_FSI = F.p1.vector()[i_f_S]
        #p_FSI = F.p_.vector()[i_f_S]
		#print "FSI pressure = ", p_FSI
		sigma_FSI = F.sigma_FSI.vector()[i_f_T]
        # structure deflection
		d_FSI  = S.d.vector()[i_s_S]
		print "structure deflection on interface = ", d_FSI

		print "fluid pressure on interface = ", p_FSI
		# fluid mesh velocity
		u_mesh_FSI = F.u_mesh.vector()[i_f_V]
		print "fluid mesh velocity on interface = ", u_mesh_FSI





        # check this indentation is correct.
        # can print mesh velocity... what about displacement? Or structure velocity. Something to compare.

        #print "fluid velocity on interface = ", u_FSI
        #print "fluid pressure on interface = ", p_FSI
		#print "about to break"
		break

	for x in F.mesh.coordinates(): x[:] += DC.dt*F.u_mesh(x)[:]
	DC.Save_Results(S, F)



# F.u1.vector() gives nodal values. Convenient to convert to standard numpy array for further processing. ie
# F.u1.vector().array()
# mesh.coordinates() returns Mxd array of coordinates. M is number of vertices and d dimension
# mesh.num_cells() returns number of cells
# mesh.num_vertices() returns number of vertices
# str(mesh) gives some mesh details.

#from numpy import where use np.where

#bc0.apply(u.vector())
#right_dofs = where(u.vector()==1)
#bc1.apply(u.vector())
#top_dofs = where(u.vector()==2)


# domains - get mesh domains, geometry
