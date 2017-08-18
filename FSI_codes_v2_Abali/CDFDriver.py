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


count = 0
i_f_d_dot = []

for i_s in range(dofs_s_V.shape[0]/2):
	for i_f in range(dofs_f_V.shape[0]):
		if (dofs_s_V[i_s*2,0] == dofs_f_V[i_f,0]) & (dofs_s_V[i_s*2,1] == dofs_f_V[i_f,1]):
			i_f_d_dot.append(i_f)
			count += count
i_f_d_dot = np.asarray(i_f_d_dot)


# locate fluid interface points that correspond to mesh
#nrows, ncols = dofs_f_V.shape
#dtype = {'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols*[dofs_f_V.dtype]}
#C = np.intersect1d(dofs_f_V.view(dtype), dofs_s_V.view(dtype), assume_unique = False)
#C = C.view(dofs_s_V.dtype).reshape(-1,ncols)


# = np.where((dofs_f_V[:,1] - dofs_s_V[:,1] == 0 ) & (dofs_f_V[:,0] - dofs_s_V[:,0] == 0 ))[0]
#ix = np.isin(dofs_f_V[:,1], dofs_s_V[:,1])
################ Iteration section ##############################
stop = 3*DC.dt
#print stop

# Sequentialy staggered iteration scheme
while DC.t < DC.T + DOLFIN_EPS:
	print ' STARTING TIME LOOP ITERATION ', DC.iter_t
	DC.t += DC.dt
	# Loop for convergence between fluid and structure
	#u_FSI = F.u1.vector()[i_f_V]
	#print "fluid velocity on interface = ", u_FSI

	for ii in range(3):
		print ''
		print ''
                # not sure what time loop iteration number is... how many to converge at given time step?
		print 'Time loop iteration number = ', DC.iter_t
		print 'Loop iteration time = ', DC.t

		# It is logical to place the fluid solver first because this causes the
		#structure to deform. However, Placing it last allows the velcities of
		# mesh, fluid and structure to be compared for the same time step.

		# Compute structural displacement and velocity
		S.Structure_Problem_Solver(DC, F)
		# Compute velocity mesh
		IO.Move_Mesh(S, F)
		# Solve fluid problem for velocity and pressure
		F.Fluid_Problem_Solver(DC, S)

	S.d0.assign(S.d) # set current time to previous for displacement
	F.u0.assign(F.u1)# set current time to previous for velocity. (pressure?)

	if DC.t == stop:
		print "Checking interface solutions at time:", stop

		# Fluid velocity, pressure and stress
		u_FSI = F.u1.vector()[i_f_V]
        #u_FSI = F.v_.vector()[i_f_V]
		# F.u1.vector().array()
		#F.u1.vector().array().shape

		# Concatenate dofs_f_V and F.u1.vector.array()
		u_map = np.concatenate((F.u1.vector().array().reshape(dofs_f_V.shape[0],1), dofs_f_V), axis = 1)

		# Calculate rate of change of deflection in order to compare to mesh velocity
		d_dot_FSI = S.d_dot.vector()[i_s_V]

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
		#print "fluid mesh velocity on interface = ", u_mesh_FSI

		#/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#engr2-1-11-24-edu:FSI_codes_v2_Abali felixnewberry$

        # check this indentation is correct.
        # can print mesh velocity... what about displacement? Or structure velocity. Something to compare.
		#print "fluid velocity on interface = ", u_FSI
		# Non zero? Not a good sign... Coordinates given by dofs_f_V[i_f_V]
        #print "fluid pressure on interface = ", p_FSI
		#print "about to break"
		print "2 norm mesh and fluid velocities on boundary:", np.linalg.norm(u_FSI-u_mesh_FSI)
		print "2 norm mesh and fluid velocites entire domain:", np.linalg.norm(F.u_mesh.vector().array() -F.u1.vector().array())
		print "2 norm mesh and structure velocites entire domain:", np.linalg.norm(F.u_mesh.vector()[i_f_d_dot] - d_dot_FSI)


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
