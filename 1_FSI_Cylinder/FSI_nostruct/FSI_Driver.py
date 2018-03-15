#  Serves as the Driver for FSI analysis of channel flow past a cylinder with an elastic bar.

# Import user-defined solvers for fluid and structure domains
from StructureSolver import StructureSolver
from FluidSolver import FluidSolver
from Cylinder import ProblemSpecific
from MeshSolver import MeshSolver

# Import dolfin and numpy
from fenics import *
import numpy as np
import pylab as plt
import scipy.io
import math

### TO DO:
# save boundary data in neater way
# save to a folder.

# consider mesh velocity in fluid variational form. Presently this has been adjusted to
# the deformed mesh, and element degree 2 function space. May affect things.

parameters["allow_extrapolation"] = True
set_log_level(ERROR)

#  Create Object with flow past cylinder expressions, mesh and subdomains
p_s = ProblemSpecific()

########### Structure Solver Setup section #################

structureElementType = "CG" # CG is equivilent to lagrange.
structureElementDegree = 1 # Standard linear lagrange element.
structureSolverMethod = "St_venant" # The only method presently available.
#StructureSolverMethod = "Linear"
structureBodyForce = Constant((0.0, 0.0))

# Set the equations for the variational solver and boundary conditions and create structure solver object
structureSolver = StructureSolver(structureElementType,
	structureElementDegree, structureSolverMethod, structureBodyForce, p_s)

########### Mesh Solver Setup section #################

meshElementType = "CG" # CG is equivilent to lagrange.
meshElementDegree = 1 # Standard linear lagrange element.

meshSolver = MeshSolver(meshElementType, meshElementDegree, p_s)

############# Fluid Solver Setup section #################

fluidElementType = "CG"		# Lagrange is the same as CG
velocityElementDegree = 2
pressureElementDegree = 1
fluidSolverMethod = "IPCS" # The only option available for the fluid solver is IPCS decomposition
fluidBodyForce = Constant((0.0,0.0))

#  Set the equations for the variational solver and boundary conditions and create fluid solver object
fluidSolver = FluidSolver(p_s.meshFluid, fluidElementType,
		velocityElementDegree, pressureElementDegree, fluidSolverMethod, fluidBodyForce, p_s, meshSolver)


########### Boundary conditions  #################

# Define boundary conditions on fluid, structure and mesh.
# Traction from fluid on the structure will be included after computing the fluid.

p_s.defineBoundaryConditions(structureSolver, fluidSolver, meshSolver, p_s)

########### Boundar Dofs  #################
# not most elegant approach

p_s.defineInterfaceDofs(structureSolver, fluidSolver, meshSolver, p_s)

#################################################################################
### Further trouble shooting things, track traction and bar velocities
#################################################################################

# Initialize tensor to record traction
# size is all nodes on FSI by x and y for fluid and structure (4) by # t steps
# number of nodes in FSI:

traction_tensor = np.zeros((p_s.dofs_s_V.shape[0],int(p_s.T/p_s.dt)+1))

# Initialize tensor to record FSI velocites.
# Size is all structure velocity nodes on FSI by 6 (mesh, fluid and structure. Split into FSI later)
#p_s.i_s_V_fsi.size
bar_vel_tensor = np.zeros((p_s.dofs_s_V.shape[0],3,int(p_s.T/p_s.dt)+1))

#### for traction comparison:
#s_normal_stresses = Function(structureSolver.scalarSpace)
#s_v=TestFunction(structureSolver.scalarSpace)

#f_normal_stresses = Function(fluidSolver.scalarSpace)
#f_v=TestFunction(fluidSolver.scalarSpace)

#################################################################################
################ Iteration section ##############################

results = np.zeros((int(p_s.T/p_s.dt)+1,4))
count = 0

#mesh displacement at same nodes as structure displacement
# all coordinates, 2 coordinates to each dof value. Want one coordinate to find two dof values. Therefore take every second.
dofs_u_disp = p_s.dofs_s_V[p_s.i_s_V_fsi[0::2]]

# matrix to store displacement. One for a given time step, one to record displacment of previous time step.
# disp_m_FSI_0 = 0*dofs_u_disp
# disp_m_FSI_1 = 0*dofs_u_disp


# load displacement matrix

# load data
# dispMatrix = np.loadtxt('artifical_disp')
# dispMatrix = np.loadtxt('artifical_disp_T16_h0_05')

# dispMatrix = np.loadtxt('artifical_disp_step_05')
# dispMatrix = np.loadtxt('artifical_disp3_T1_h_05_dt05')
# dispMatrix = np.loadtxt('artifical_disp_med2_h_005_stept2')

# try load just the t vector and the disp vector separately.



# dispMatrix = np.loadtxt('artifical_disp_step_001')
# dispMatrix = 0*dispMatrix

# requre an additional column of zeros to set initial... maybe can avoid this.
# size is 294 by 5001
# to select just one column
# dispMatrix[:,n]
#self.u0.vector()[:] = u_temp
# structureSolver.u12, structureSolver.v12 = structureSolver.U.split(deepcopy = True)
# structureSolver.meshInterfaceDisplacement = interpolate(structureSolver.u12,structureSolver.interfaceDispVectorSpace)

# meshInterfaceDisplacement = Function(structureSolver.interfaceDispVectorSpace)
# disp_temp0 = np.ascontiguousarray(dispMatrix[:,0], dtype = np.float64)
# disp_temp1 = np.ascontiguousarray(dispMatrix[:,0], dtype = np.float64)

# dispParabola = np.loadtxt('artifical_disp_med2_h_005_parabola')
# dispTime = np.loadtxt('artifical_disp_med2_h_005_time')

struct_coords = np.loadtxt('struct_coords')
x_min = min(struct_coords[:,0])
x_max = max(struct_coords[:,0])
x_length = x_max - x_min;

height = 0.005;

period = 1;
dispParabola = height/x_length**2.0*(struct_coords[:,0]-x_min)**2;
dispParabola[0::2] = 0 
# set y values to 0.

# dispTime = sin((2.0*pi)/(period)*p_s.t)


disp_temp0 = np.ascontiguousarray(0*dispParabola, dtype = np.float64)
disp_temp1 = np.ascontiguousarray(0*dispParabola, dtype = np.float64)


# Sequentialy staggered iteration scheme
while p_s.t + p_s.dt<= p_s.T: # + DOLFIN_EPS:
	print 'STARTING TIME LOOP ITERATION ', count

	p_s.t += p_s.dt

	# artificially set structure displacement
	# disp_temp1 = np.ascontiguousarray(dispMatrix[:,count], dtype = np.float64)
	# disp_temp1 = np.ascontiguousarray(dispTime[count]*dispParabola, dtype = np.float64)
	disp_temp1 = np.ascontiguousarray(sin((2.0*pi)/(period)*p_s.t)*dispParabola, dtype = np.float64)

	# structureSolver.meshInterfaceDisplacement.vector()[:] = disp_temp

	structureSolver.meshInterfaceDisplacement.vector()[:] = disp_temp1

	for ii in range(3):
		print ''
		print ''
		print 'Time loop iteration number = ', ii
		print 'Loop iteration time = ', p_s.t

		# Loop for convergence between fluid and structure, updated mesh velocity
        # changes fluid solution which updates traction on structure which updates
		# velocity on mesh


		# Boundary vector function fluidInterfaceVelocity2 breaks down for some coordinates at a given time step approx 108
		# It does not break when only mesh is done.
		# write analogous mesh  here to compare to.

		# Interpolate or project degree 1 to degree 2.

		fluidSolver.fluidInterfaceVelocityTest.vector()[:] = meshSolver.meshInterfaceVelocity.vector().get_local()

		# This is the function used in BC which seems to break after approx 80 time steps
		fluidSolver.fluidInterfaceVelocity2Test = project(fluidSolver.fluidInterfaceVelocityTest, fluidSolver.interfaceFluidVectorSpace2, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )

		# Solve fluid problem for velocity and pressure
		fluidSolver.solveFluidProblem(p_s, meshSolver.meshInterfaceVelocity, meshSolver.meshDisplacement, meshSolver.meshVelocity)

		# fluid velocity should match structure velocity which should match mesh velocity.
		# fluid velocity on FSI at nodes that also appear on structure...
		u_f_FSI = fluidSolver.u1.vector()[p_s.i_f_V_fsi_com]
		u_m_FSI = meshSolver.meshVelocity.vector()[p_s.i_m_V_fsi_com]

		# Structure displacement and velocity on FSI
		u11, v11 = structureSolver.U.split(deepcopy = True)
		u10, v10 = structureSolver.U0.split(deepcopy = True)

		# u_s_FSI = u11.vector()[p_s.i_s_V_fsi]
		u_s_FSI = disp_temp1

		# v_s_FSI = v11.vector()[p_s.i_s_V_fsi]

# interfaceStructureDispCoords = structureSolver.interfaceStructureDispCoords
# these coordinates
# p_s.dofs_m_V[p_s.i_m_V_fsi_com] - interfaceStructureDispCoords

		v_s_FSI = structureSolver.u_vel.vector()[p_s.i_s_V_fsi]
		# doesn't match nodes so ignore for now...

		v_s_FSI = (1.0/p_s.dt)*(disp_temp1 - disp_temp0)


		# 294 displacements. use coordinates to find mesh displacements at same pionts
		dispCoords = p_s.dofs_s_V[p_s.i_s_V_fsi]
		# u_s_FSI is half lenth of dispCoords (two coordinates for each velocity)

		# does not match structure velocity well. one is midpoint, other is not.
		# could compare to structure disp1 - structure disp2...

		# u_m_FSI nodes 140 to 150 seem suspect.

		# p_s.i_m_V_fsi_com[140:150]
		# p_s.dofs_m_V[p_s.i_m_V_fsi_com[140:150]]
		# u_m_FSI[141:148]
		# u_f_FSI[141:148]
		# p_s.dofs_m_V[p_s.i_m_V_fsi_com[141:148]]

		# fluid boundary set from expression which is based off this vector:
		#fluidSolver.fluidInterfaceVelocity.vector.get_local() has 0 values in many places, as expected in x direction... it should be non zero
		#where mesh is non zero for x = 0.6 points.
		# the vector on which expression is based appears to have the correct mesh values.
		# therefore, assuming the expression is applied correction, the mistaken values are because the bc is not applied to the end of the bar.
		# This could potentially explain the slight discrepency in traction integrals too.
		# look at medium mesh facet files...
		# have to recall mesh conversion code...

		# interfaceFluidCoords = fluidSolver.interfaceFluidVectorSpace.tabulate_dof_coordinates().reshape((fluidSolver.interfaceFluidVectorSpace.dim(),-1))

		# 141 to 148 again...
		interfaceFluidCoords2 = fluidSolver.interfaceFluidVectorSpace.tabulate_dof_coordinates().reshape((fluidSolver.interfaceFluidVectorSpace.dim(),-1))

		interfaceMeshCoords2 = meshSolver.interfaceMeshVectorSpace.tabulate_dof_coordinates().reshape((meshSolver.interfaceMeshVectorSpace.dim(),-1))

		print "2 norm mesh and fluid velocities :", np.linalg.norm(u_m_FSI - u_f_FSI)/np.linalg.norm(u_m_FSI)
		print "2 norm mesh and structure velocities :", np.linalg.norm(u_m_FSI - v_s_FSI)/np.linalg.norm(u_m_FSI)
		# print "2 norm fluid and structure velocities :", np.linalg.norm(u_f_FSI - v_s_FSI)/np.linalg.norm(u_f_FSI)

		# # hmm, traction integral on boundary does differ a smidge
		# print "integral of fluid traction:", fluidSolver.integral_0
		# print "integral of structure traction:", fluidSolver.integral_1

		# displacement is not even close... how does that work? Ah, updating on every step...

		# disp_m_FSI = 0*u_s_FSI
		# for i_disp in range(u_s_FSI.size/2):
		# 	# x
		# 	disp_m_FSI[2*i_disp] = meshSolver.u1(dispCoords[2*i_disp])[0]
		# 	# y
		# 	disp_m_FSI[2*i_disp+1] = meshSolver.u1(dispCoords[2*i_disp])[1]
		#
		# print "2 norm mesh and structure displacement :", np.linalg.norm(u_s_FSI - disp_m_FSI)/np.linalg.norm(u_s_FSI)

		# Solve structure problem for displacement and velocity
		# structureSolver.structureProblemSolver(p_s, fluidSolver)

		# t_f_FSI = fluidSolver.traction.vector()[p_s.i_m_V_fsi_com]

		# in structure solver traction is applied with an expression.
		# structureSolver.fluidInterfaceTractionExpression
		# tractionCoords = p_s.dofs_m_V[p_s.i_m_V_fsi_com]

		# t_s_FSI = 0*t_f_FSI

		# for i_traction in range(tractionCoords.size/4):
			# step through every second coordinate. (size of coords is 2x length) append x and then y values
			# x value
			# t_s_FSI[2*i_traction] = structureSolver.fluidInterfaceTractionExpression(tractionCoords[2*i_traction])[0]
			# y value
			# t_s_FSI[2*i_traction+1] = structureSolver.fluidInterfaceTractionExpression(tractionCoords[2*i_traction])[1]


        # structureSolver.fluidInterfaceTractionExpression(p_s.dofs_m_V[p_s.i_m_V_fsi_com])
		# Traction on interface
		# one should be negative of the other.
		# print "2 norm fluid and structure traction :", np.linalg.norm(t_f_FSI - t_s_FSI)/np.linalg.norm(t_f_FSI)
		# has a couple of large values...


		# compare fluid, structure and mesh velocity.
		# mesh and structure displacement
		# fluid and structure traction.

		# Solve mesh problem
		meshSolver.meshProblemSolver(disp_temp1, fluidSolver, p_s)

		disp_m_FSI = 0*u_s_FSI
		for i_disp in range(u_s_FSI.size/2):
			# x
			disp_m_FSI[2*i_disp] = meshSolver.u1(dispCoords[2*i_disp])[0]
			# y
			disp_m_FSI[2*i_disp+1] = meshSolver.u1(dispCoords[2*i_disp])[1]

		print "2 norm mesh and structure displacement :", np.linalg.norm(u_s_FSI - disp_m_FSI)/np.linalg.norm(u_s_FSI)

		# Update fluid mesh coordinates
		for i_coords in range(len(fluidSolver.mesh.coordinates())):
			fluidSolver.mesh.coordinates()[i_coords] = meshSolver.mesh.coordinates()[i_coords] + meshSolver.u1(meshSolver.mesh.coordinates()[i_coords])

		for i_coords in range(len(p_s.interfaceFluid.coordinates())):
			p_s.interfaceFluid.coordinates()[i_coords] = p_s.interfaceMesh.coordinates()[i_coords] + meshSolver.meshInterfaceDisplacement(p_s.interfaceMesh.coordinates()[i_coords])

		# Midpoint:
		# for i_coords in range(len(fluidSolver.mesh.coordinates())):
			# fluidSolver.mesh.coordinates()[i_coords] = meshSolver.mesh.coordinates()[i_coords] + meshSolver.meshDisplacement(meshSolver.mesh.coordinates()[i_coords])

		# for i_coords in range(len(p_s.interfaceFluid.coordinates())):
			# p_s.interfaceFluid.coordinates()[i_coords] = p_s.interfaceMesh.coordinates()[i_coords] + meshSolver.meshInterfaceDisplacement(p_s.interfaceMesh.coordinates()[i_coords])



	# update FSI displacement
	# disp_m_FSI_0 = disp_m_FSI_1

	# Update values

	# Update structure solver
	#structureSolver.u, structureSolver.v = structureSolver.U.split()

	#split(deepcopy = True)
	# propogate displacements and velocities
	#structureSolver.u1, structureSolver.v1 = structureSolver.U.split()
	# This one line should be sufficient update for structure. u0 and v0 are read from U0
	# structureSolver.U0.assign(structureSolver.U)
	#structureSolver.u0.assign(structureSolver.u1)
	#structureSolver.v0.assign(structureSolver.v1)
	#structureSolver.u0, structureSolver.v0 = structureSolver.U0.split(deepcopy = True)
	#structureSolver.u0, structureSolver.v0 = split(structureSolver.U0)

	fluidSolver.u0.assign(fluidSolver.u1)	# set current time to previous for velocity.
	fluidSolver.p0.assign(fluidSolver.p1)

	# disp_temp0.assign(disp_temp1)
	disp_temp0 = disp_temp1

	meshSolver.u0.assign(meshSolver.u1)
	#print "2 norm mesh and structure displacements :", np.linalg.norm(u_f_FSI_t - u_s_FSI)

	## Compute drag and lift
	#drag, lift = p_s.compute_forces(p_s.mesh_f,p_s.nu_f, fluidSolver.u1, fluidSolver.p1, fluidSolver.ds)

	#print "cylinder drag and lift = ", drag, lift

	p_s.saveResults(structureSolver, fluidSolver, meshSolver)

	#Lift and drag forces acting on the whole submerged body
	# (cylinder and structure together)
	drag, lift = p_s.computeForces(p_s.meshFluid,p_s.mu_f, fluidSolver.u1, fluidSolver.p1, fluidSolver.ds)

	#print " x and y displacements at point A:", A_disp

	# print "ux, uy, drag and lift on cylinder and bar:", [u11(0.6,0.2)[0], u11(0.6,0.2)[1],drag, lift]

	# results[count,:] = [A_disp_s.l_s + p_s.D/2. + p_s.lp[0], A_disp[1], drag, lift]

	# compute and save traction
	# Compute x and y traction in structure

	# Structure
	# s_fx=(1/FacetArea(structureSolver.mesh))*s_v*structureSolver.traction[0]*structureSolver.ds(20)
	# s_fy=(1/FacetArea(structureSolver.mesh))*s_v*structureSolver.traction[1]*structureSolver.ds(20)

	# s_Traction_x=assemble(s_fx,tensor=s_normal_stresses.vector())[p_s.i_s_S_fsi]
	# s_Traction_y=assemble(s_fy,tensor=s_normal_stresses.vector())[p_s.i_s_S_fsi]

	# Fluid
	# n_f = FacetNormal(fluidSolver.mesh)
	# T_f = dot(fluidSolver.sigma_FSI, n_f)

	# f_fx=(1/FacetArea(fluidSolver.mesh))*f_v*T_f[0]*fluidSolver.ds(20)
	# f_fy=(1/FacetArea(fluidSolver.mesh))*f_v*T_f[1]*fluidSolver.ds(20)

	# f_Traction_x=assemble(f_fx,tensor=f_normal_stresses.vector())[i_f_S_fsi]
	# f_Traction_y=assemble(f_fy,tensor=f_normal_stresses.vector())[i_f_S_fsi]

	# size i_f_S_fsi.size, 4, dt
	# traction_tensor[:,:,count] = np.column_stack((f_Traction_x,f_Traction_y,s_Traction_x,s_Traction_y))
	# traction_tensor[:,count] = t_s_FSI

	# bar_vel_tensor[:, :, count] = np.column_stack((v_s_FSI, u_m_FSI, u_f_FSI))

	##u11, v11 = structureSolver.U.split(deepcopy = True)
	# results[count,:] = [u11(0.6,0.2)[0], u11(0.6,0.2)[1], drag, lift]

	##if count%500 == 0:
		##scipy.io.savemat('cylinder_fsi_course_1.mat', mdict={'results':results})

	# # move mesh
	# problematic because although x is one thing in fluid mesh, it is another in reference mesh
	# likely that meshVelocity(x) extrapolates.

	# for x in fluidSolver.mesh.coordinates(): x[:] += p_s.dt*meshSolver.meshVelocity(x)[:]

	# update fluid mesh and fluid FSI boundary mesh with new displacements
	# relative to reference mesh.



	# fluid velocity coords are transferred on boundary mesh. Check coordinates of this, compare to reference.
	# p_s.interfaceMesh.coordinates()[70:84]
	# p_s.interfaceFluid.coordinates()
	# mesh_diff = p_s.interfaceFluid.coordinates()-p_s.interfaceMesh.coordinates()
	# print out reference coordinates.
		# meshSolver.mesh.coordinates()[i_coords]


	# update mesh displacement accordingly
	# update mesh displacement using mesh velocity.
	# first, set previos mesh displacement to present
	# meshSolver.u0.vector()[:] = meshSolver.meshDisplacement1.vector().get_local()

	# this should only occur when time step changes!!! not within time step.. ie displacement remains at 0 until mesh is moved.
	# Update present
	# meshSolver.meshDisplacement1.vector()[:] += p_s.dt*meshSolver.meshVelocity.vector().get_local()

	# take average over time step for use in piola transform of fluid stress
	# may have to initialize meshDisplacement.
	# meshSolver.meshDisplacement.vector()[:] = 0.5*(meshSolver.meshDisplacement0.vector().get_local()+ meshSolver.meshDisplacement1.vector().get_local())

	# check that mesh displacement is close to structure displacement


	print count

	count += 1

	if count%500 == 0:
		# Save FSI values and validation metrics

		pwd_restart = './Restart_FSI/'
		pwd_results = './Results_FSI/'

		## compute and save nodal values
		nodal_values_u = fluidSolver.u1.vector().get_local()
		np.savetxt(pwd_restart+'nodal_u_1', nodal_values_u)
		nodal_values_p = fluidSolver.p1.vector().get_local()
		np.savetxt(pwd_restart+'nodal_p_1', nodal_values_p)
		nodal_values_m = meshSolver.meshVelocity.vector().get_local()
		np.savetxt(pwd_restart+'nodal_m', nodal_values_m)
		nodal_values_disp = meshSolver.u1.vector().get_local()
		np.savetxt(pwd_restart+'nodal_m_disp1', nodal_values_disp)
		nodal_values_disp = meshSolver.u0.vector().get_local()
		np.savetxt(pwd_restart+'nodal_m_disp0', nodal_values_disp)

		# try saving just the U
		#u, v = structureSolver.U.split()
		nodal_values_U = structureSolver.U.vector().get_local()
		np.savetxt(pwd_restart+'nodal_Ustr_1', nodal_values_U)

		scipy.io.savemat(pwd_results+'cylinder_fsi_course.mat', mdict={'results':results})

		scipy.io.savemat(pwd_results+'traction_mesh_med.mat', mdict={'traction_tensor':traction_tensor})
		scipy.io.savemat(pwd_results+'coordinates_fsi_traction_vector.mat', mdict={'tractionCoords':tractionCoords})

		scipy.io.savemat(pwd_results+'bar_vel_tensor_mesh_med.mat', mdict={'bar_vel_tensor':bar_vel_tensor})

		coordinates_fsi_Vector = p_s.dofs_s_V[p_s.i_s_V_fsi]
		scipy.io.savemat(pwd_results+'coordinates_fsi_mesh_med.mat', mdict={'coordinates_fsi_Vector':coordinates_fsi_Vector})


	print("t =", p_s.t)


# Save FSI values and validation metrics

pwd_restart = './Restart_FSI/'
pwd_results = './Results_FSI/'


## compute and save nodal values
nodal_values_u = fluidSolver.u1.vector().get_local()
np.savetxt(pwd_restart+'nodal_u_1', nodal_values_u)
nodal_values_p = fluidSolver.p1.vector().get_local()
np.savetxt(pwd_restart+'nodal_p_1', nodal_values_p)
nodal_values_m = meshSolver.meshVelocity.vector().get_local()
np.savetxt(pwd_restart+'nodal_m', nodal_values_m)
nodal_values_disp = meshSolver.u1.vector().get_local()
np.savetxt(pwd_restart+'nodal_m_disp1', nodal_values_disp)
nodal_values_disp = meshSolver.u0.vector().get_local()
np.savetxt(pwd_restart+'nodal_m_disp0', nodal_values_disp)


# np.savetxt('struct_coords',structureSolver.interfaceStructureDispCoords)


np.savetxt(pwd_restart+'nodal_m_disp0', nodal_values_disp)

# try saving just the U
#u, v = structureSolver.U.split()
nodal_values_U = structureSolver.U.vector().get_local()
np.savetxt(pwd_restart+'nodal_Ustr_1', nodal_values_U)


scipy.io.savemat(pwd_results+'cylinder_fsi_stiff.mat', mdict={'results':results})

scipy.io.savemat(pwd_results+'traction_mesh_med.mat', mdict={'traction_tensor':traction_tensor})
scipy.io.savemat(pwd_results+'coordinates_fsi_traction_vector.mat', mdict={'tractionCoords':tractionCoords})

scipy.io.savemat(pwd_results+'bar_vel_tensor_mesh_med.mat', mdict={'bar_vel_tensor':bar_vel_tensor})

coordinates_fsi_Vector = p_s.dofs_s_V[p_s.i_s_V_fsi]
scipy.io.savemat(pwd_results+'coordinates_fsi_mesh_med.mat', mdict={'coordinates_fsi_Vector':coordinates_fsi_Vector})

# save coordinates and vector values for structure displacement.
interfaceStructureDispCoords = structureSolver.interfaceStructureDispCoords
scipy.io.savemat(pwd_results+'disp_coordinates.mat', mdict={'interfaceStructureDispCoords':interfaceStructureDispCoords})

# disp vector to see what it looks like:
nodal_disp = structureSolver.meshInterfaceDisplacement.vector().get_local()
scipy.io.savemat(pwd_results+'nodal_disp.mat', mdict={'nodal_disp':nodal_disp})


### Some trouble shooting script

# pwd = './Results_Cylinder_trouble_shooting/'


# # save fluid before movement
#
# # save mesh before movement
# file_u_m_before = File(pwd + 'u_f.pvd')
# file_u_m_before << (meshSolver.u1)
#
# # move mesh
#
# # first set displacement
# disp_temp1 = np.ascontiguousarray(dispMatrix[:,5], dtype = np.float64)
# # structureSolver.meshInterfaceDisplacement.vector()[:] = disp_temp
#
# structureSolver.meshInterfaceDisplacement.vector()[:] = disp_temp1
# meshSolver.meshProblemSolver(disp_temp1, fluidSolver, p_s)
#
# # update mesh coordiantes
# for i_coords in range(len(fluidSolver.mesh.coordinates())):
# 	fluidSolver.mesh.coordinates()[i_coords] = meshSolver.mesh.coordinates()[i_coords] + meshSolver.u1(meshSolver.mesh.coordinates()[i_coords])
#
# for i_coords in range(len(p_s.interfaceFluid.coordinates())):
# 	p_s.interfaceFluid.coordinates()[i_coords] = p_s.interfaceMesh.coordinates()[i_coords] + meshSolver.meshInterfaceDisplacement(p_s.interfaceMesh.coordinates()[i_coords])
#
# # save fluid after movement - perhaps this will not see updated mesh
# file_u_f_after = File(pwd + 'u_f1.pvd')
# file_u_f_after << (fluidSolver.u1)
#
# # save mesh displacement
# file_u_m = File(pwd + 'u_m.pvd')
# file_u_m << (meshSolver.u1)
# # save mesh velocity
# file_v_m = File(pwd + 'v_m.pvd')
# file_v_m << (meshSolver.meshVelocity)
#
#
# fluidSolver.meshVelocityCurrent1.vector()[:] = meshSolver.meshVelocity.vector().get_local()
#
# # Interpolate or project degree 1 to degree 2.
# meshLocal = project(fluidSolver.meshVelocityCurrent1, fluidSolver.vectorSpace, solver_type = "mumps", \
# 	form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )
#
# # Compute fluid velocity - mesh velocity and save This
#
# u_augmented = fluidSolver.u1
# u_augmented.vector()[:] = fluidSolver.u1.vector().get_local() - meshLocal.vector().get_local()
#
# # save mesh local
# # save mesh velocity
# file_v_m_local = File(pwd + 'v_m_local.pvd')
# file_v_m_local << (meshLocal)
#
# file_u_augmented = File(pwd + 'u_augmented.pvd')
# file_u_augmented << (u_augmented)

# Peclet number calculation


# # DG scalar space for Peclet number
# pecletSpace  = FunctionSpace(fluidSolver.mesh, "DG", 1)
# pecletNumber = Function(pecletSpace2)
#
# # map cells to dofs
# dofmap = pecletSpace2.dofmap()
# dof_to_cell = [cell for cell in range(fluidSolver.mesh.num_cells())
#                for dof in dofmap.cell_dofs(cell)]
#
# # step through dofs, identify appropriate cell, compute Pe
# dofs_peclet = pecletSpace2.tabulate_dof_coordinates().reshape((pecletSpace2.dim(),-1))
#
# for i_coords in range(len(dofs_peclet)):
# 	# dofs_peclet[i_coords]
# 	u_mag = sqrt(fluidSolver.u1(dofs_peclet[i_coords])[0]**2 + fluidSolver.u1(dofs_peclet[i_coords])[1]**2)
#
# 	cell_index = dof_to_cell[i_coords]
# 	cell_1 = Cell(fluidSolver.mesh, cell_index)
# 	cell_size = cell_1.circumradius()
#
# 	# Pe_dof = u_mag*cell_size/p_s.nu_f
# 	pecletNumber.vector()[i_coords] = u_mag*cell_size/p_s.nu_f
#
# # Save Peclet number
# pwd = './Results_Cylinder_trouble_shooting/'
# file_Pe = File(pwd + 'Pe.pvd')
# file_Pe << (pecletNumber)



# Could I step through coordinates of cells instead?

# step through coordinates
# dir(fluidSolver.mesh)

	# store Pe for all cells... save to file?
	# could save a Peclet number for each vertex... then take average. Yeah.
# calc norm of fluid velocity

# Load reference mesh on which mesh eqn is solved.


### find coordinates for artificial displacement.
meshStruct2 = Mesh('meshFiles/cylinderbar_fine2.xml');
facetsMeshStruct2 = MeshFunction("size_t", meshStruct2, "meshFiles/cylinderbar_med2_facet_region.xml")

# restrict mesh to exterior boundaries only:
boundaryStructure2 = BoundaryMesh(meshStruct2, 'exterior')
interfaceStructure2 = SubMesh(boundaryStructure2, p_s.fSIInterface)

interfaceDispVectorSpace2 = VectorFunctionSpace(interfaceStructure2, "CG", 1)

interfaceStructureDispCoords2 = interfaceDispVectorSpace2.tabulate_dof_coordinates().reshape((interfaceDispVectorSpace2.dim(),-1))

# save coordinates and vector values for structure displacement.
# interfaceStructureDispCoords = structureSolver.interfaceStructureDispCoords
scipy.io.savemat(pwd_results+'disp_coordinates3.mat', mdict={'interfaceStructureDispCoords2':interfaceStructureDispCoords2})
