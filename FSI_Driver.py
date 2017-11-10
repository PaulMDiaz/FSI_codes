#  Serves as the Driver for FSI analysis of Cavity
#  Driven Flow with a bottom non-rigid wall
#
#  Created: 26 February 2014
#  Modified: 26 July 2016

# Import user-defined solvers for fluid and structure domains
from StructureSolver import Structure_Solver
from FluidSolver import Fluid_Solver
from Cylinder import problem_specific

from MeshSolver import *

# Import dolfin and numpy
from fenics import *
import numpy as np
import pylab as plt
import scipy.io
import math

parameters["allow_extrapolation"] = True
set_log_level(ERROR)

#  Create Object with Driven Cavity expressions, mesh and subdomains
p_s = problem_specific()

Mesh_Solver = Mesh_Solver()

############# Fluid Solver Setup section #################
#  Set the Fluid Element Type
#FluidElementType = "Lagrange"
FluidElementType = 'CG'
# Lagrange is the same as CG

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
F = Fluid_Solver(p_s.mesh_f, FluidElementType,
		VelocityElementDegree, PressureElementDegree, FluidSolverMethod, FluidBodyForce)

########### Structure Solver Setup section #################

#  Set the Structure Element Type
#StructureElementType = "CG" # CG is equivilent to lagrange.
StructureElementType = "Lagrange" # CG is equivilent to lagrange.

#  Set the Structure Element Degree
StructureElementDegree = 1 # Standard linear lagrange element.

# Set the solver used for the structure problem
# The options are "Linear" or "NeoHookean" for the structure solver
#StructureSolverMethod = "NeoHookean"
StructureSolverMethod = "St_venant"
#StructureSolverMethod = "Linear"
# Body forces on the structure
StructureBodyForce = Constant((0.0, 0.0))

# Define boundary conditions on subdomains (fluid and structure) but traction from fluid on the structure which
# will be included after computing the fluid.


# Set the equations for the variational solver and boundary conditions and creat structure solver object
S = Structure_Solver(p_s.mesh_s, StructureElementType,
	StructureElementDegree, StructureSolverMethod, StructureBodyForce)

p_s.Define_Boundary_Conditions(S, F)


### for printing the interaction

# Fluid  array of coordinates

# S scalar, V vector, T Tensor
dofs_f_S = F.S_space.tabulate_dof_coordinates().reshape((F.S_space.dim(),-1))
dofs_f_V = F.V_space.tabulate_dof_coordinates().reshape((F.V_space.dim(),-1))
dofs_f_T = F.T_space.tabulate_dof_coordinates().reshape((F.T_space.dim(),-1))

# Structure array of coordinates
dofs_s_S = S.S_space.tabulate_dof_coordinates().reshape((S.S_space.dim(),-1))
dofs_s_V = S.V_space.tabulate_dof_coordinates().reshape((S.V_space.dim(),-1))
dofs_s_T = S.T_space.tabulate_dof_coordinates().reshape((S.T_space.dim(),-1))

dofs_s_V2 = S.V2_space.tabulate_dof_coordinates().reshape((S.V2_space.dim(),-1))

#Extract dof indices for values on boundary.
# y = 0.5 if mesh is not deformed.
#i_f_S = np.where((dofs_f_S[:,1] == p_s.h))[0] #  & (x <= 0.5)

# inlet
i_f_V_in = np.where((dofs_f_V[:,0] == 0.0))[0] #  & (x <= 0.5)

# monitor Point
i_s_A = np.where((dofs_s_V[:,0] == p_s.x_bar + p_s.l) & (dofs_s_V[:,1] == 0.19 + p_s.h/2) )[0]

i_s_A2 = np.where((dofs_s_V2[:,0] == p_s.x_bar + p_s.l) & (dofs_s_V2[:,1] == 0.19 + p_s.h/2) )[0]
#i_f_T = np.where((dofs_f_T[:,1] == p_s.h))[0] #  & (x <= 0.5)

i_f_Q_1 = np.where((dofs_f_S[:,0] == 0.2-0.05))[0]
i_f_Q_2 = np.where((dofs_f_S[:,0] == 0.2+0.05))[0]
#i_s_S = np.where((dofs_s_S[:,1] == p_s.h))[0] #  & (x <= 0.5)
#i_s_V = np.where((dofs_s_V[:,1] == p_s.h))[0] #  & (x <= 0.5)
#_s_T = np.where((dofs_s_T[:,1] == p_s.h))[0] #  & (x <= 0.5)
#
# Extract dofdofs_f_S indices for fluid velocity on top plate
#i_f_V_top = np.where((dofs_f_V[:,1] == p_s.H))[0]
#
# identify coordinates on structure boundary that match coordinates on fluid boundary.
# Find indices in fluid vector space.

#i_f_scomp = []

#for i_s in range(dofs_s_V[i_s_V].shape[0]/2):
#	for i_f in range(dofs_f_V[i_f_V].shape[0]/2):
#		if dofs_s_V[i_s_V[2*i_s]][0] == dofs_f_V[i_f_V[2*i_f]][0]:
#			i_f_scomp.append(i_f_V[2*i_f])
#			i_f_scomp.append(i_f_V[2*i_f+1])

################ Iteration section ##############################

# Initialize array to store A displacements, lift and drag. Should converge to steady state.

results = np.zeros((int(p_s.T/p_s.dt)+1,4))
count = 0
# Sequentialy staggered iteration scheme
while p_s.t < p_s.T + DOLFIN_EPS:
	print 'STARTING TIME LOOP ITERATION ', p_s.iter_t
	p_s.t += p_s.dt


	for ii in range(3): #change to 3 to run properly.
		print ''
		print ''
                # not sure what time loop iteration number is... how many to converge at given time step?
		print 'Time loop iteration number = ', p_s.iter_t
		print 'Loop iteration time = ', p_s.t

		# It is logical to place the fluid solver first because this causes the
		#structure to deform. However, placing it last allows the velcities of
		# mesh, fluid and structure to be compared for the same time step.

		# Loop for convergence between fluid and structure
		#u_FSI = F.u1.vector()[i_f_V]
		#print "fluid velocity on interface = ", u_FSI

		# Solve fluid problem for velocity and pressure
		F.Fluid_Problem_Solver(p_s, S)

		# Fluid velocity should match structure and mesh of previous time step.
		#u_inlet = F.u_.vector()[i_f_V_in]
		#print "fluid velocity at inlet = ", u_inlet

		#d_FSI  = S.d.vector()[i_s_S]
		#print "structure deflection on interface prior to structure solve = ", d_FSI
		# Compute structural displacement and velocity

		#pressure_f = F.p1.vector()[i_f_S]
		#print 'pressure_FSI = ', pressure_f

		# Print stress on interface from fluid
		#sigma_FSI = F.sigma_FSI.vector()[i_f_T]
		#print 'Fluid sigma_FSI = ', sigma_FSI

		# Structure velocity of previous time step should match fluid velcity. Maybe.
		#if p_s.t >= 2*p_s.dt or ii >= 2:
		#	d_FSI_1  = S.d.vector()[i_s_V]
		#	# structure velocity
		#	d_dot_FSI = S.d_dot.vector()[i_s_V]

		# Solve structure problem for displacement and displancement rate of change
		S.Structure_Problem_Solver(p_s, F)

		### Save Sigma_FSI for use in other codes
		#nodal_F_sigmaFSI = F.sigma_FSI.vector().array()
		#np.savetxt('nodal_F_sigmaFSI', nodal_F_sigmaFSI)
#
	#	nodal_S_sigmaFSI = S.sigma_FSI.vector().array()
		#np.savetxt('nodal_S_sigmaFSI', nodal_S_sigmaFSI)



		#if p_s.t >= 2*p_s.dt or ii >= 2:
			#u_M_FSI_1 = F.u_mesh.vector()[i_f_V]
			#u_M_FSI_2 = F.u_mesh.vector()[i_f_scomp]

		# Compute mesh velocity
		Mesh_Solver.Move_Mesh(S, F)

		# Mesh displacement on interface.

		# The structure and mesh velocites should match when after mesh solver is computed
		#if p_s.t >= 2*p_s.dt or ii >= 2:
		#	print "2 norm mesh and fluid velocities :", np.linalg.norm(u_M_FSI_1 - u_f)
		#	print "2 norm mesh and structure velocities :", np.linalg.norm(u_M_FSI_2 - d_dot_FSI)
			#print "2 norm fluid and structure velocities :", np.linalg.norm(u_f - d_dot_FSI)

		#d_FSI_1 = S.d.vector()[i_s_V]
		#print "structure displacement on interface", d_FSI_1

	# Compute and store monitor point
	#print S.d_.vector().array()
	# creates two new functions u and v...
	u1, v1 = S.U.split(deepcopy = True)
	A_disp = u1.vector()[i_s_A]

	# Lift and drag forces acting on the whole submerged body
	# (cylinder and structure together)
	drag, lift = p_s.compute_forces(F.mesh, F.p_, F.ds)

	#p_diff = F.p_.vector()[i_f_Q_1]-F.p_.vector()[i_f_Q_2]

	print " x and y displacements at point A:", A_disp

	print "drag and lift on cylinder and bar:", [drag, lift]

	results[count,:] = [A_disp[0], A_disp[1], drag, lift]

	count += 1
	# update structure solver
	# propogate displacements and velocities
	S.U0.assign(S.U)

	# extract displacements and velocities for results...
	#u, v = S.U.split()



	#S.d_n.assign(S.d_) 	# set current time to previous for displacement

	F.u_n.assign(F.u_)	# set current time to previous for velocity.
	#S.d00_s.assign(S.d0)

	p_s.Save_Results(S,F)
	#for x in F.mesh.coordinates(): x[:] += p_s.dt*F.u_mesh(x)[:]
	#p_s.Save_Results(S, F)

# Record y velocity along mid y line. Record x velocity for x  = 0
#u_v = np.column_stack((np.linspace(0.0, p_s.W, 101), np.zeros(101)))
#u_u = np.column_stack((np.linspace(p_s.h, p_s.H, 101), np.zeros(101)))

#for i_x in range(len(u_v)):
	#u_v[i_x,1] = F.u1(u_v[i_x,0],(p_s.H -p_s.h)/2+p_s.h)[1]
	#u_u[i_x,1] = F.u1(p_s.W/2,u_u[i_x,0])[0]

#scipy.io.savemat('u_v.mat', mdict={'u_v':u_v})
#scipy.io.savemat('u_u.mat', mdict={'u_u':u_u})

scipy.io.savemat('results_N128_t0001.mat', mdict={'results':results})

## compute and save nodal values
nodal_values_u = F.u_.vector().array()
np.savetxt('nodal_u', nodal_values_u)
nodal_values_p = F.p_.vector().array()
np.savetxt('nodal_p', nodal_values_p)
nodal_values_m = F.u_mesh.vector().array()
np.savetxt('nodal_m', nodal_values_m)
# try saving just the U
#u1, v1 = S.U.split(deepcopy = True)
#nodal_values_u1 = u1.vector().array()
#np.savetxt('nodal_u1', nodal_values_u1)
#nodal_values_v1 = v1.vector().array()
#np.savetxt('nodal_v1', nodal_values_v1)
nodal_values_U = S.U.vector().array()
np.savetxt('nodal_U', nodal_values_U)
