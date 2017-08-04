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

################ Iteration section ##############################

# Sequentialy staggered iteration scheme
while DC.t < DC.T + DOLFIN_EPS:
	print ' STARTING TIME LOOP ITERATION ', DC.iter_t
	DC.t += DC.dt
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


                # Code check, exit while loop after first time step.
                # Print fluid stress, sturcture and mesh displacement along interface. 
        if DC.t == DC.dt:
                print 'Checking interface solutions at time:', DC.dt
                # interface_values boundary is FSI() marked as F.fsi.mark(F.facets, 3)
                #print "Fluid stuff: %d" %(F.sigma_FSI)
                #sys.exit() # stops the script.
                #print "Fluid stuff: %d" %( )
                #u_array = F.u1
                u_nodal_array = F.u1.vector().array()
                print u_nodal_array[42]
                #print u_nodal_array.shape
                #print F.u1.vector().shape
                # Get vertices sitting on boundary
                #d2v = dof_to_vertex_map(F.S_space)
                #vertices_on_boundary = d2v
                #print "fluid velocity on interface = " , F.u1.vector()[F.v_.vector() == 3].array()
                #print "fluid velocity on interface = " , F.u1.vector()[F.v_.vector() == ]

                #print F.mesh.num_vertices()
                #print F.mesh.num_cells()
                #print len(u_nodal_array)

                # number of coordinates is less then u soltutions
                
                #coor = F.mesh.coordinates()
                #if F.mesh.num_vertices() == len(u_nodal_array):
                #        for i in range(F.mesh.num_vertices()):
                #                print 'u1(%g, %g) = %g' % (coor[i][0], coor[i][1], u_nodal_array[i])
                F.fsi
                break
        
               
                
	S.d0.assign(S.d) # set current time to previous for displacement
	F.u0.assign(F.u1)# set current time to previous for velocity. (pressure?)
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
