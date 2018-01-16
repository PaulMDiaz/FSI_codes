#  Serves as the Driver for FSI analysis of channel flow past a cylinder with an elastic bar.

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

#  Create Object with flow past cylinder expressions, mesh and subdomains
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
# The only option available for the fluid solver is IPCS decomposition
FluidSolverMethod = "IPCS"

# Set the body forces on the fluid
FluidBodyForce = Constant((0.0,0.0))

#  Set the equations for the variational solver and boundary conditions and create fluid solver object
fluid_solver = Fluid_Solver(p_s.mesh_f, FluidElementType,
		VelocityElementDegree, PressureElementDegree, FluidSolverMethod, FluidBodyForce, p_s)

########### Structure Solver Setup section #################

#  Set the Structure Element Type
#StructureElementType = "CG" # CG is equivilent to lagrange.
StructureElementType = "CG" # CG is equivilent to lagrange.

#  Set the Structure Element Degree
StructureElementDegree = 1 # Standard linear lagrange element.

# Set the solver used for the structure problem
# The options are "Linear" or "NeoHookean" for the structure solver
#StructureSolverMethod = "NeoHookean"
StructureSolverMethod = "St_venant"
#StructureSolverMethod = "Linear"
# Body forces on the structure
StructureBodyForce = Constant((0.0, 0.0))

# Set the equations for the variational solver and boundary conditions and create structure solver object
S = Structure_Solver(p_s.mesh_s, StructureElementType,
	StructureElementDegree, StructureSolverMethod, StructureBodyForce, p_s, fluid_solver)

# Define boundary conditions on subdomains (fluid and structure). Traction from fluid on the structure
# will be included after computing the fluid.

p_s.Define_Boundary_Conditions(S, fluid_solver)

### for printing the interaction. Some details here are really problem specific.
# Should eventually be commented out once code is sound.

# Key
# S scalar, V vector, T Tensor. f fluid, s structure
dofs_f_S = fluid_solver.S_space.tabulate_dof_coordinates().reshape((fluid_solver.S_space.dim(),-1))
dofs_f_V = fluid_solver.V_space.tabulate_dof_coordinates().reshape((fluid_solver.V_space.dim(),-1))
dofs_f_T = fluid_solver.T_space.tabulate_dof_coordinates().reshape((fluid_solver.T_space.dim(),-1))

# Structure array of coordinates
dofs_s_S = S.S_space.tabulate_dof_coordinates().reshape((S.S_space.dim(),-1))
dofs_s_V = S.V_space.tabulate_dof_coordinates().reshape((S.V_space.dim(),-1))
dofs_s_T = S.T_space.tabulate_dof_coordinates().reshape((S.T_space.dim(),-1))

# Structure with mixed element function space
#dofs_s_V2 = S.V_space.tabulate_dof_coordinates().reshape((S.V_space.dim(),-1))

#Extract dof indices for pints of of interest.
# monitor Point - end of bar
i_s_A = np.where((dofs_s_V[:,0] == p_s.x_bar + p_s.l) & (dofs_s_V[:,1] == 0.19 + p_s.h/2) )[0]
#i_s_A2 = np.where((dofs_s_V2[:,0] == p_s.x_bar + p_s.l) & (dofs_s_V2[:,1] == 0.19 + p_s.h/2) )[0]

# find indices for FSI
# Both fluid and structure, scalar and vector indices
# Break interface into three regions, top, bottom and end.
# Scalar
# top
i_s_S_fsi1 = np.where((dofs_s_S[:,1] == 0.21) & (dofs_s_S[:,0] >= p_s.x_bar) & (dofs_s_S[:,0] <= 0.6))[0]
i_f_S_fsi1 = np.where((dofs_f_S[:,1] == 0.21) & (dofs_f_S[:,0] >= p_s.x_bar) & (dofs_f_S[:,0] <= 0.6))[0]
# bottom
i_s_S_fsi2 = np.where((dofs_s_S[:,1] == 0.19) & (dofs_s_S[:,0] >= p_s.x_bar) & (dofs_s_S[:,0] <= 0.6))[0]
i_f_S_fsi2 = np.where((dofs_f_S[:,1] == 0.19) & (dofs_f_S[:,0] >= p_s.x_bar) & (dofs_f_S[:,0] <= 0.6))[0]
# end
i_s_S_fsi3 = np.where((dofs_s_S[:,0] == 0.6) & (dofs_s_S[:,1] > 0.19) & (dofs_s_S[:,1] < 0.21))[0]
i_f_S_fsi3 = np.where((dofs_f_S[:,0] == 0.6) & (dofs_f_S[:,1] > 0.19) & (dofs_f_S[:,1] < 0.21))[0]

# Reorder and concatenate such that indices number nodes left to right along the bottom,
# up the end and then right to left along the top.

i_s_S_fsi1 = i_s_S_fsi1[dofs_s_S[i_s_S_fsi1][:,0].argsort()][::-1]
i_f_S_fsi1 = i_f_S_fsi1[dofs_f_S[i_f_S_fsi1][:,0].argsort()][::-1]
i_s_S_fsi2 = i_s_S_fsi2[dofs_s_S[i_s_S_fsi2][:,0].argsort()]
i_f_S_fsi2 = i_f_S_fsi2[dofs_f_S[i_f_S_fsi2][:,0].argsort()]
i_s_S_fsi3 = i_s_S_fsi3[dofs_s_S[i_s_S_fsi3][:,1].argsort()]
i_f_S_fsi3 = i_f_S_fsi3[dofs_f_S[i_f_S_fsi3][:,1].argsort()]

# Concatenate bottom then end then top. 1 and 2 are in a good order.
i_s_S_fsi = np.concatenate((i_s_S_fsi2,i_s_S_fsi3,i_s_S_fsi1),axis = 0)
i_f_S_fsi = np.concatenate((i_f_S_fsi2,i_f_S_fsi3,i_f_S_fsi1),axis = 0)

# Vector
# top
i_s_V_fsi1 = np.where((dofs_s_V[:,1] == 0.21) & (dofs_s_V[:,0] >= p_s.x_bar) & (dofs_s_V[:,0] <= 0.6))[0]
i_f_V_fsi1 = np.where((dofs_f_V[:,1] == 0.21) & (dofs_f_V[:,0] >= p_s.x_bar) & (dofs_f_V[:,0] <= 0.6))[0]
# bottom
i_s_V_fsi2 = np.where((dofs_s_V[:,1] == 0.19) & (dofs_s_V[:,0] >= p_s.x_bar) & (dofs_s_V[:,0] <= 0.6))[0]
i_f_V_fsi2 = np.where((dofs_f_V[:,1] == 0.19) & (dofs_f_V[:,0] >= p_s.x_bar) & (dofs_f_V[:,0] <= 0.6))[0]
# end
i_s_V_fsi3 = np.where((dofs_s_V[:,0] == 0.6) & (dofs_s_V[:,1] >= 0.19 +DOLFIN_EPS) & (dofs_s_V[:,1] <= 0.21 - DOLFIN_EPS))[0]
i_f_V_fsi3 = np.where((dofs_f_V[:,0] == 0.6) & (dofs_f_V[:,1] >= 0.19+DOLFIN_EPS) & (dofs_f_V[:,1] <= 0.21 -DOLFIN_EPS))[0]

# Seperate into x and y components so that reordering works fine.

# Reorder and concatenate such that indices number nodes left to right along the bottom,
# up the end and then right to left along the top.
#i_s_V_fsi1 = i_s_V_fsi1[dofs_s_V[i_s_V_fsi1][:,0].argsort()]#[::-1]
#i_f_V_fsi1 = i_f_V_fsi1[dofs_f_V[i_f_V_fsi1][:,0].argsort()]#[::-1]
# break into pairs
#i_s_V_fsi1_pairs = zip(i_s_V_fsi1[::2], i_s_V_fsi1[1::2])
#i_f_V_fsi1_pairs = zip(i_f_V_fsi1[::2], i_f_V_fsi1[1::2])
# reverse pair order
#i_f_V_fsi1_pairs = i_f_V_fsi1_pairs[::-1]
#i_s_V_fsi1_pairs = i_s_V_fsi1_pairs[::-1]
# put into list
#i_f_V_fsi1 = list(sum(i_f_V_fsi1_pairs, ()))
#i_s_V_fsi1 = list(sum(i_s_V_fsi1_pairs, ()))

#i_s_V_fsi2 = i_s_V_fsi2[dofs_s_V[i_s_V_fsi2][:,0].argsort()]
#i_f_V_fsi2 = i_f_V_fsi2[dofs_f_V[i_f_V_fsi2][:,0].argsort()]
#i_s_V_fsi3 = i_s_V_fsi3[dofs_s_V[i_s_V_fsi3][:,1].argsort()]
#i_f_V_fsi3 = i_f_V_fsi3[dofs_f_V[i_f_V_fsi3][:,1].argsort()]


# Concatenate
i_s_V_fsi = np.concatenate((i_s_V_fsi2,i_s_V_fsi3,i_s_V_fsi1),axis = 0)
i_f_V_fsi = np.concatenate((i_f_V_fsi2,i_f_V_fsi3,i_f_V_fsi1),axis = 0)

# This is a complete set of indices. Find ones that match between fluid and structure:
# This might not be the best approach. For instance, interpolating the structure solution.

# Fluid has more nodes because for lagrange elements of order greater than 2 there are
# nodes that do not correspond to vertices.

# overlapping nodes between structure and fluid
i_f_V_fsi_com = []

for i_s in range(dofs_s_V[i_s_V_fsi].shape[0]/2):
	for i_f in range(dofs_f_V[i_f_V_fsi].shape[0]/2):
		if dofs_s_V[i_s_V_fsi[2*i_s]][0] == dofs_f_V[i_f_V_fsi[2*i_f]][0] and dofs_s_V[i_s_V_fsi[2*i_s]][1] == dofs_f_V[i_f_V_fsi[2*i_f]][1]:
			i_f_V_fsi_com.append(i_f_V_fsi[2*i_f])
			i_f_V_fsi_com.append(i_f_V_fsi[2*i_f+1])


################ Iteration section ##############################

# Initialize tensor to record traction
# size is all nodes on FSI by x and y for fluid and structure (4) by # t steps
traction_tensor = np.zeros((i_f_S_fsi.size,4,int(p_s.T/p_s.dt)+1))

# for traction comparison:
s_normal_stresses = Function(S.S_space)
s_v=TestFunction(S.S_space)

f_normal_stresses = Function(fluid_solver.S_space)
f_v=TestFunction(fluid_solver.S_space)

################ Iteration section ##############################

results = np.zeros((int(p_s.T/p_s.dt)+1,4))
count = 0

#mesh displacement at same nodes as structure displacement
# all coordinates, 2 coordinates to each dof value. Want one coordinate to find two dof values. Therefore take every second.
dofs_u_disp = dofs_s_V[i_s_V_fsi[0::2]]

# matrix to store displacement. One for a given time step, one to record displacment of previous time step.
disp_m_FSI_0 = 0*dofs_u_disp
disp_m_FSI_1 = 0*dofs_u_disp

# Sequentialy staggered iteration scheme
while p_s.t <= p_s.T: # + DOLFIN_EPS:
	print 'STARTING TIME LOOP ITERATION ', count

	p_s.t += p_s.dt

	for ii in range(3):
		print ''
		print ''
		print 'Time loop iteration number = ', ii
		print 'Loop iteration time = ', p_s.t

		# Loop for convergence between fluid and structure.
		# Solve fluid problem for velocity and pressure
		fluid_solver.solve_fluid_problem(p_s)

		# fluid velocity should match structure velocity which should match mesh velocity.
		# fluid velocity on FSI at nodes that also appear on structure...
		u_f_FSI = fluid_solver.u1.vector()[i_f_V_fsi_com]

		# Structure displacement and velocity on FSI

		u_s_FSI = S.u1.vector()[i_s_V_fsi]
		v_s_FSI = S.v1.vector()[i_s_V_fsi]

		# mesh velocity on FSI
		u_m_FSI = fluid_solver.u_mesh.vector()[i_f_V_fsi_com]

		#print "2 norm mesh and fluid velocities :", np.linalg.norm(u_m_FSI/alpha_1 - u_f_FSI/alpha_1)
		#print "2 norm mesh and structure velocities :", np.linalg.norm(u_m_FSI/alpha_2 - v_s_FSI/alpha_2)
		#print "2 norm fluid and structure velocities :", np.linalg.norm(u_f_FSI/alpha_3 - v_s_FSI/alpha_3)
		#print "2 norm mesh and structure deflections :", np.linalg.norm(u_s_FSI/alpha_4 - disp_m_FSI_0.flatten()/alpha_4)

		print "2 norm mesh and fluid velocities :", np.linalg.norm(u_m_FSI - u_f_FSI)
		print "2 norm mesh and structure velocities :", np.linalg.norm(u_m_FSI - v_s_FSI)
		print "2 norm fluid and structure velocities :", np.linalg.norm(u_f_FSI - v_s_FSI)
		print "2 norm mesh and structure deflections :", np.linalg.norm(u_s_FSI - disp_m_FSI_0.flatten())

		#np.linalg.norm(dofs_s_V[i_s_V_fsi] - dofs_f_V[i_f_V_fsi_com])

		# Solve structure problem for displacement and velocity
		S.Structure_Problem_Solver(p_s, fluid_solver)

		# Solve mesh problem
		Mesh_Solver.Move_Mesh(S, fluid_solver)

	# update FSI displacement
	disp_m_FSI_1 = disp_m_FSI_0

	# Update values

	# Update structure solver
	#S.u, S.v = S.U.split()

	#split(deepcopy = True)
	# propogate displacements and velocities
	#S.u1, S.v1 = S.U.split()
	#S.U0.assign(S.U)
	#S.u0.assign(S.u1)
	#S.v0.assign(S.v1)
	#S.u0, S.v0 = S.U0.split(deepcopy = True)
	#S.u0, S.v0 = split(S.U0)

	# Propagate the displacements, velocities and accelerations
	S.u0.assign(S.u1)
	S.v0.assign(S.v1)
	S.a0.assign(S.a1)


	fluid_solver.u0.assign(fluid_solver.u1)	# set current time to previous for velocity.
	fluid_solver.p0.assign(fluid_solver.p1)

	#print "2 norm mesh and structure displacements :", np.linalg.norm(u_f_FSI_t - u_s_FSI)

	## Compute drag and lift
	#drag, lift = p_s.compute_forces(p_s.mesh_f,p_s.nu_f, fluid_solver.u1, fluid_solver.p1, fluid_solver.ds)

	#print "cylinder drag and lift = ", drag, lift

	p_s.Save_Results(S,fluid_solver)

	#A_disp = u1.vector()[i_s_A]

	#Lift and drag forces acting on the whole submerged body
	# (cylinder and structure together)
	drag, lift = p_s.compute_forces(p_s.mesh_f,p_s.nu_f, fluid_solver.u1, fluid_solver.p1, fluid_solver.ds)

	#print " x and y displacements at point A:", A_disp

	print "ux, uy, drag and lift on cylinder and bar:", [S.u1(0.6,0.2)[0], S.u1(0.6,0.2)[1],drag, lift]

	#results[count,:] = [A_disp_s.l_s + p_s.D/2. + p_s.lp[0], A_disp[1], drag, lift]

	# compute and save traction
	# Compute x and y traction in structure

	# Structure
	s_fx=(1/FacetArea(S.mesh))*s_v*S.T_hat[0]*S.ds(2)
	s_fy=(1/FacetArea(S.mesh))*s_v*S.T_hat[1]*S.ds(2)

	s_Traction_x=assemble(s_fx,tensor=s_normal_stresses.vector())[i_s_S_fsi]
	s_Traction_y=assemble(s_fy,tensor=s_normal_stresses.vector())[i_s_S_fsi]

	# Fluid
	n_f = FacetNormal(fluid_solver.mesh)
	T_f = dot(fluid_solver.sigma_FSI, n_f)

	f_fx=(1/FacetArea(fluid_solver.mesh))*f_v*T_f[0]*fluid_solver.ds(3)
	f_fy=(1/FacetArea(fluid_solver.mesh))*f_v*T_f[1]*fluid_solver.ds(3)

	f_Traction_x=assemble(f_fx,tensor=f_normal_stresses.vector())[i_f_S_fsi]
	f_Traction_y=assemble(f_fy,tensor=f_normal_stresses.vector())[i_f_S_fsi]

	# size i_f_S_fsi.size, 4, dt
	traction_tensor[:,:,count-1] = np.column_stack((f_Traction_x,f_Traction_y,s_Traction_x,s_Traction_y))

	results[count,:] = [u1(0.6,0.2)[0], u1(0.6,0.2)[1], drag, lift]


	# move mesh
	for x in fluid_solver.mesh.coordinates(): x[:] += p_s.dt*fluid_solver.u_mesh(x)[:]

	count += 1

	#u_fsi_t = fluid_solver.u1.vector()[i_f_V_fsi_t]

	# things go bad when u_mesh is used as fluid BC...

	#u_inlet = fluid_solver.u1.vector()[i_f_V_in]
	#print u_inlet
	#for x in fluid_solver.mesh.coordinates(): x[:] += p_s.dt*fluid_solver.u_mesh(x)[:]
	#p_s.Save_Results(S, fluid_solver)

# Record y velocity along mid y line. Record x velocity for x  = 0
#u_v = np.column_stack((np.linspace(0.0, p_s.W, 101), np.zeros(101)))
#u_u = np.column_stack((np.linspace(p_s.h, p_s.H, 101), np.zeros(101)))

#for i_x in range(len(u_v)):
	#u_v[i_x,1] = fluid_solver.u1(u_v[i_x,0],(p_s.H -p_s.h)/2+p_s.h)[1]
	#u_u[i_x,1] = fluid_solver.u1(p_s.W/2,u_u[i_x,0])[0]

#scipy.io.savemat('u_v.mat', mdict={'u_v':u_v})
#scipy.io.savemat('u_u.mat', mdict={'u_u':u_u})

## compute and save nodal values
#nodal_values_u = fluid_solver.u1.vector().array()
#np.savetxt('nodal_u_32', nodal_values_u)
#nodal_values_p = fluid_solver.p1.vector().array()
#np.savetxt('nodal_p_32', nodal_values_p)
#nodal_values_m = fluid_solver.u_mesh.vector().array()
#np.savetxt('nodal_m_32', nodal_values_m)3
# try saving just the U
#u, v = S.U.split()
#nodal_values_U = S.U.vector().array()
#np.savetxt('nodal_U_32', nodal_values_U)

#scipy.io.savemat('results_64.mat', mdict={'results':results})

scipy.io.savemat('traction_tensor.mat', mdict={'traction_tensor':traction_tensor})

scipy.io.savemat('u_m_FSI.mat', mdict={'u_m_FSI':u_m_FSI})
scipy.io.savemat('v_s_FSI.mat', mdict={'v_s_FSI':v_s_FSI})
scipy.io.savemat('u_f_FSI.mat', mdict={'u_f_FSI':u_f_FSI})
