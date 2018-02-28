
from dolfin import *
import numpy as np
import pylab as plt
from mshr import *

class ProblemSpecific:

	def __init__(self):

        ############## INPUT DATA PARAMETERS ###################
	    # Physical parameters
		self.nu_f = 0.001	# Fluid kinematic viscosity
		self.nu_s = 0.4	# Structure Poisson coefficient
		self.mu_s = 0.5e6 # structure first lame constant (very small)
		self.mu_s = 0.5e11 # structure first lame constant (very stiff, to test)

		self.rho_f = 1000.0	# Fluid density (incorporated in the fluid corrected pressure as p_corr = p/rho)
		self.rho_s = 1000.0
		self.mu_f = self.rho_f*self.nu_f # dynamic visosity

		self.E_s = self.mu_s*2*(1+self.nu_s)	# structure elastic modulus
		self.lambda_s = self.nu_s*self.E_s/((1+self.nu_s)*(1-2*self.nu_s)) # 2nd lame constant

		self.E_s = Constant(self.E_s)
		self.mu_s = Constant(self.mu_s)
		self.nu_s = Constant(self.nu_s)
		self.rho_s = Constant(self.rho_s)

		# mean inlet velocity
		self.U_mean = 0.2

		# Numerical parameters
		self.dt = 0.001 # Time step
		self.T = 8.000  #  Set final time for iteration

		# Geometric parameters
		self.L = 2.5 	#Channel length
		self.H = 0.41	#Channel height
		self.l = 0.35 	#Bar length
		self.h = 0.02	#Bar height

		# Set up a variable for time
		self.t_vec = np.linspace(0.0,self.T,self.T/self.dt+1)
		self.t = self.t_vec[0]
		self.iter_t = 1

		# load meshes
		# med, fine or course
		self.meshStructure = Mesh('meshFiles/cylinderbar_med.xml');
		self.facetsStructure = MeshFunction("size_t", self.meshStructure, "meshFiles/cylinderbar_med_facet_region.xml")

		# self.meshStructure = Mesh('cylinderbar_course.xml');
		# self.facetsStructure = MeshFunction("size_t", self.meshStructure, "cylinderbar_course_facet_region.xml")

		#self.sub_domains_structure = MeshFunction("size_t", self.meshStructure, "cylinderbar_physical_region.xml")
		self.meshFluid = Mesh('meshFiles/cylinderfluid_med.xml');
		self.facetsFluid = MeshFunction("size_t", self.meshFluid, "meshFiles/cylinderfluid_med_facet_region.xml")
		# self.meshFluid = Mesh('cylinderfluid_course.xml');
		# self.facetsFluid = MeshFunction("size_t", self.meshFluid, "cylinderfluid_course_facet_region.xml")

		# Load reference mesh on which mesh eqn is solved.
		self.meshRef = Mesh('meshFiles/cylinderfluid_med.xml');
		self.facetsRef = MeshFunction("size_t", self.meshRef, "meshFiles/cylinderfluid_med_facet_region.xml")

		# numbering system:
 		# subdomans: 1 fluid, 1 Structure (seperated by resepective subsomains..)
		# boundaries 15 walls, 16 inlet, 17 outlet, 18 cylinder, 29 fixed bar end, 20 Fsi

		# Define boundary meshes
		# Set up top, end and bottom. Inelegant number for greater than, but works
		# class FSIInterface(SubDomain):
		#     def inside(self, x, on_boundary):
		#         return x[0] > 0.24898979 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS \
		# 		and x[1] > 0.19 - DOLFIN_EPS and x[1] < 0.19 + DOLFIN_EPS \
		# 		or x[0] > 0.24898979 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS \
		# 		and x[1] < 0.21 +DOLFIN_EPS  and x[1] > 0.21 - DOLFIN_EPS \
		# 		or x[0] > 0.6 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS \
		# 		and x[1] < 0.21 + DOLFIN_EPS and x[1] > 0.19 - DOLFIN_EPS \

		class FSIInterface(SubDomain):
		    def inside(self, x, on_boundary):
		        return x[0] > 0.24 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS \
				and x[1] > 0.19 - DOLFIN_EPS and x[1] < 0.19 + DOLFIN_EPS \
				or x[0] > 0.24 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS \
				and x[1] < 0.21 + DOLFIN_EPS  and x[1] > 0.21 - DOLFIN_EPS \
				or x[0] > 0.6 - DOLFIN_EPS and x[0] < 0.6 + DOLFIN_EPS \
				and x[1] < 0.21 + DOLFIN_EPS and x[1] > 0.19 - DOLFIN_EPS \

		self.fSIInterface = FSIInterface()

		#################
		# mesh
		#################
		boundariesMesh = MeshFunction("size_t", self.meshRef, self.meshRef.topology().dim() - 1)
		boundariesMesh.set_all(0)
		self.fSIInterface.mark(boundariesMesh, 1)

		# restrict mesh to exterior boundaries only:
		boundaryMesh = BoundaryMesh(self.meshRef, 'exterior')
		self.interfaceMesh = SubMesh(boundaryMesh, self.fSIInterface)

		#################
		# fluid
		#################
		self.boundariesFluid = MeshFunction("size_t", self.meshFluid, self.meshFluid.topology().dim() - 1)
		self.boundariesFluid.set_all(0)
		self.fSIInterface.mark(self.boundariesFluid, 1)

		# restrict mesh to exterior boundaries only:
		boundaryFluid = BoundaryMesh(self.meshFluid, 'exterior')
		self.interfaceFluid = SubMesh(boundaryFluid, self.fSIInterface)

		#################
		# structure
		#################

		# Not necessary. Traction is calculated over entire domain and set as an expression.

		#self.mesh = generate_mesh(self.domain, self.N)
		# self.A_point = Point(0.6, 0.2)

		L = self.L
		H = self.H
		h = self.h
		l = self.l

		U_mean = self.U_mean

		# Variables to generate files
		pwd = './Results_Cylinder_FSI_med_trouble/'
		self.file_u_s = File(pwd + 'u_s.pvd')
		self.file_v_s = File(pwd + 'v_s.pvd')
		self.file_v_f = File(pwd + 'v_f.pvd')
		self.file_p_f = File(pwd + 'p_f.pvd')
		self.file_v_m = File(pwd + 'v_m.pvd')

	def computeForces(self,Mesh,mu,u,p, ds):
		self.mesh = Mesh
		self.mu = mu
		# uses both deformed mesh and normal of deformed mesh.

		# Face normals
		n = FacetNormal(self.mesh)
		# Stress tensor
		# Traction
		sigma = self.mu*(grad(u)+grad(u).T) -p*Identity(2)
		T = dot(sigma,n)

		self.drag = -T[0]*ds(18)-T[0]*ds(20)
		self.lift = -T[1]*ds(18)-T[1]*ds(20)

		self.drag1 = assemble(self.drag); self.lift1 = assemble(self.lift);
		return self.drag1, self.lift1

 # boundaries already defined in geometry file...
	def defineBoundaryConditions(self, structureSolver, fluidSolver, meshSolver, p_s):

		L = self.L
		H = self.H
		h = self.h
		l = self.l

		U_mean = self.U_mean

		# numbering system:
 		# subdomans: 1 fluid, 1 Structure (seperated by resepective subsomains..)
		# boundaries 15 walls, 16 inlet, 17 outlet, 18 cylinder, 29 bar end, 20 Fsi

		############# Initialize structure boundary #############
		structureSolver.cells = CellFunction("size_t", structureSolver.mesh)
		structureSolver.ds = Measure('ds', domain = structureSolver.mesh, subdomain_data = p_s.facetsStructure)
		structureSolver.dx = Measure('dx', domain = structureSolver.mesh, subdomain_data = structureSolver.cells)
		structureSolver.n = FacetNormal(structureSolver.mesh)

		#  BCs for the left side (no displacement)
		LeftBC = DirichletBC(structureSolver.mixedSpace.sub(0), Constant((0.0, 0.0)), p_s.facetsStructure, 19)
		#  Set up the boundary conditions
		structureSolver.bcs = [LeftBC]
		# FSI boundary is addressed in variational form with Neuman condition

		############ Initialize fluid boundary ###################
		# Create mesh function and mark all 0
		fluidSolver.cells = CellFunction("size_t", fluidSolver.mesh)
		fluidSolver.ds = Measure('ds', domain = fluidSolver.mesh, subdomain_data = p_s.facetsFluid)
		fluidSolver.dx = Measure('dx', domain = fluidSolver.mesh, subdomain_data = fluidSolver.cells)
		fluidSolver.n = FacetNormal(fluidSolver.mesh)

		inlet_profile = ('1.5*U_mean*x[1]*(H-x[1])/pow(H/2, 2)', '0')
		bcuInlet = DirichletBC(fluidSolver.vectorSpace, Expression(inlet_profile, H = H, U_mean = U_mean, degree=2), p_s.facetsFluid, 16)
		bcuWalls = DirichletBC(fluidSolver.vectorSpace, Constant((0, 0)), p_s.facetsFluid, 15)
		bcuCylinder = DirichletBC(fluidSolver.vectorSpace, Constant((0, 0)), p_s.facetsFluid, 18)

		# define expression for FSI boundary
		class FluidInterfaceVelocityExpression(Expression):
		    def eval(self, values, x):
		        try:
		            values[:] = fluidSolver.fluidInterfaceVelocity(x)
		        except:
		            values[:] = 0
		    def value_shape(self):
		        return (2,)

		fluidInterfaceVelocityExpression = FluidInterfaceVelocityExpression(degree=2)

		bcuFSI = DirichletBC(fluidSolver.vectorSpace, fluidInterfaceVelocityExpression, p_s.facetsFluid, 20)

		# Pressure
		bcpOutlet = DirichletBC(fluidSolver.scalarSpace, Constant(0), p_s.facetsFluid, 17)

		# Set up the boundary conditions
		fluidSolver.bcu = [bcuInlet, bcuWalls, bcuCylinder, bcuFSI]
		fluidSolver.bcp = [bcpOutlet]

		############ Initialize mesh boundary ###################
		meshSolver.cells = CellFunction("size_t", meshSolver.mesh)
		meshSolver.ds = Measure('ds', domain = meshSolver.mesh, subdomain_data = p_s.facetsRef)
		meshSolver.dx = Measure('dx', domain = meshSolver.mesh, subdomain_data = meshSolver.cells)
		meshSolver.n = FacetNormal(meshSolver.mesh)

		# Extract velocity to use in mesh solver
		u12, v12 = structureSolver.U.split()
		u12.set_allow_extrapolation(True)
		v12.set_allow_extrapolation(True)

		bcMeshInlet = DirichletBC(meshSolver.vectorSpace, Constant((0,0)), p_s.facetsRef, 16)
		bcMeshWalls = DirichletBC(meshSolver.vectorSpace, Constant((0,0)), p_s.facetsRef, 15)
		bcMeshCylinder = DirichletBC(meshSolver.vectorSpace, Constant((0,0)), p_s.facetsRef, 18)
		bcMeshOutlet = DirichletBC(meshSolver.vectorSpace, Constant((0,0)), p_s.facetsRef, 17)

		# should be displacement.
		bcMeshFSI = DirichletBC(meshSolver.vectorSpace, u12, p_s.facetsRef, 20)

		meshSolver.bcMesh = [bcMeshInlet, bcMeshWalls, bcMeshCylinder, bcMeshOutlet, bcMeshFSI]

	def defineInterfaceDofs(self, structureSolver, fluidSolver, meshSolver, p_s):

		# Mesh solution
		dofsMesh = meshSolver.vectorSpace.tabulate_dof_coordinates().reshape((meshSolver.vectorSpace.dim(),-1))
		# Mesh interface
		dofsMeshInterface = meshSolver.interfaceMeshVectorSpace.tabulate_dof_coordinates().reshape((meshSolver.interfaceMeshVectorSpace.dim(),-1))
		# Fluid interface
		dofsFluidInterface = fluidSolver.interfaceFluidVectorSpace.tabulate_dof_coordinates().reshape((fluidSolver.interfaceFluidVectorSpace.dim(),-1))
		# Fluid solution
		dofsFluid = fluidSolver.vectorSpace.tabulate_dof_coordinates().reshape((fluidSolver.vectorSpace.dim(),-1))

		# mesh values on boundary:
		indices_mesh = np.where((dofsMesh[:,1] > 0.19-DOLFIN_EPS) & (dofsMesh[:,1] < 0.21 +DOLFIN_EPS) & (dofsMesh[:,0] > 0.24898979-DOLFIN_EPS) & (dofsMesh[:,0] < 0.6+DOLFIN_EPS))[0]
		# a different order between interface and mesh. should be okay with interpolate.
		dofsMesh[indices_mesh]

		# Key
		# S scalar, V vector, T Tensor. f fluid, s structure
		dofs_f_S = fluidSolver.scalarSpace.tabulate_dof_coordinates().reshape((fluidSolver.scalarSpace.dim(),-1))
		dofs_f_V = fluidSolver.vectorSpace.tabulate_dof_coordinates().reshape((fluidSolver.vectorSpace.dim(),-1))
		dofs_f_T = fluidSolver.tensorSpace.tabulate_dof_coordinates().reshape((fluidSolver.tensorSpace.dim(),-1))

		# S scalar, V vector, T Tensor. f fluid, s structure
		self.dofs_m_V = meshSolver.vectorSpace.tabulate_dof_coordinates().reshape((meshSolver.vectorSpace.dim(),-1))

		# Structure array of coordinates
		dofs_s_S = structureSolver.scalarSpace.tabulate_dof_coordinates().reshape((structureSolver.scalarSpace.dim(),-1))
		self.dofs_s_V = structureSolver.vectorSpace.tabulate_dof_coordinates().reshape((structureSolver.vectorSpace.dim(),-1))
		dofs_s_T = structureSolver.tensorSpace.tabulate_dof_coordinates().reshape((structureSolver.tensorSpace.dim(),-1))

		# Structure with mixed element function space
		dofs_s_V2 = structureSolver.mixedSpace.tabulate_dof_coordinates().reshape((structureSolver.mixedSpace.dim(),-1))

		#Extract dof indices for pints of of interest.
		# monitor Point - end of bar
		i_s_A = np.where((self.dofs_s_V[:,0] == 0.6) & (self.dofs_s_V[:,1] == 0.19 + p_s.h/2) )[0]
		i_s_A2 = np.where((dofs_s_V2[:,0] == 0.6 + p_s.l) & (dofs_s_V2[:,1] == 0.19 + p_s.h/2) )[0]

		# find indices for FSI
		# need more generic
		i_s_S_fsi1 = np.where((dofs_s_S[:,1] == 0.21) & (dofs_s_S[:,0] >= 0.24) & (dofs_s_S[:,0] <= 0.6))[0]
		i_s_S_fsi1 = np.where((dofs_s_S[:,1] <= 0.21 + DOLFIN_EPS) & (dofs_s_S[:,1] >= 0.21 - DOLFIN_EPS))[0]

		i_f_S_fsi1 = np.where((dofs_f_S[:,1] == 0.21) & (dofs_f_S[:,0] >= 0.24) & (dofs_f_S[:,0] <= 0.6))[0]
		i_f_S_fsi1 = np.where((dofs_f_S[:,1] <= 0.21 + DOLFIN_EPS) & (dofs_f_S[:,1] >= 0.21 - DOLFIN_EPS))[0]

		# bottom
		i_s_S_fsi2 = np.where((dofs_s_S[:,1] == 0.19) & (dofs_s_S[:,0] >= 0.24) & (dofs_s_S[:,0] <= 0.6))[0]
		i_s_S_fsi2 = np.where((dofs_s_S[:,1] <= 0.19 + DOLFIN_EPS) & (dofs_s_S[:,1] >= 0.19 - DOLFIN_EPS))[0]

		i_f_S_fsi2 = np.where((dofs_f_S[:,1] == 0.19) & (dofs_f_S[:,0] >= 0.24) & (dofs_f_S[:,0] <= 0.6))[0]
		i_f_S_fsi2 = np.where((dofs_f_S[:,1] <= 0.19 + DOLFIN_EPS) & (dofs_f_S[:,1] >= 0.19 - DOLFIN_EPS))[0]

		# end
		i_s_S_fsi3 = np.where((dofs_s_S[:,0] == 0.6) & (dofs_s_S[:,1] > 0.19) & (dofs_s_S[:,1] < 0.21))[0]
		i_s_S_fsi3 = np.where((dofs_s_S[:,0] <= 0.6 + DOLFIN_EPS) & (dofs_s_S[:,0] >= 0.6 - DOLFIN_EPS)& (dofs_s_S[:,1] > 0.19) &(dofs_s_S[:,1] < 0.21 - DOLFIN_EPS))[0]

		i_f_S_fsi3 = np.where((dofs_f_S[:,0] == 0.6) & (dofs_f_S[:,1] > 0.19) & (dofs_f_S[:,1] < 0.21))[0]
		i_f_S_fsi3 = np.where((dofs_f_S[:,0] <= 0.6 + DOLFIN_EPS) & (dofs_f_S[:,0] >= 0.6 - DOLFIN_EPS)& (dofs_f_S[:,1] > 0.19) &(dofs_f_S[:,1] < 0.21 - DOLFIN_EPS))[0]

		# Reorder and concatenate such that indices number nodes left to right along the bottom,
		# up the end and then right to left along the top.

		i_s_S_fsi1 = i_s_S_fsi1[dofs_s_S[i_s_S_fsi1][:,0].argsort()][::-1]
		i_f_S_fsi1 = i_f_S_fsi1[dofs_f_S[i_f_S_fsi1][:,0].argsort()][::-1]
		i_s_S_fsi2 = i_s_S_fsi2[dofs_s_S[i_s_S_fsi2][:,0].argsort()]
		i_f_S_fsi2 = i_f_S_fsi2[dofs_f_S[i_f_S_fsi2][:,0].argsort()]
		i_s_S_fsi3 = i_s_S_fsi3[dofs_s_S[i_s_S_fsi3][:,1].argsort()]
		i_f_S_fsi3 = i_f_S_fsi3[dofs_f_S[i_f_S_fsi3][:,1].argsort()]

		# Concatenate bottom then end then top. 1 and 2 are in a good order.
		self.i_s_S_fsi = np.concatenate((i_s_S_fsi2,i_s_S_fsi3,i_s_S_fsi1),axis = 0)
		i_f_S_fsi = np.concatenate((i_f_S_fsi2,i_f_S_fsi3,i_f_S_fsi1),axis = 0)

		# Vector
		# top
		i_s_V_fsi1 = np.where((self.dofs_s_V[:,1] == 0.21) & (self.dofs_s_V[:,0] >= 0.24) & (self.dofs_s_V[:,0] <= 0.6))[0]
		i_f_V_fsi1 = np.where((dofs_f_V[:,1] == 0.21) & (dofs_f_V[:,0] >= 0.24) & (dofs_f_V[:,0] <= 0.6))[0]
		i_m_V_fsi1 = np.where((self.dofs_m_V[:,1] == 0.21) & (self.dofs_m_V[:,0] >= 0.24) & (self.dofs_m_V[:,0] <= 0.6))[0]

		# bottom
		i_s_V_fsi2 = np.where((self.dofs_s_V[:,1] == 0.19) & (self.dofs_s_V[:,0] >= 0.24) & (self.dofs_s_V[:,0] <= 0.6))[0]
		i_f_V_fsi2 = np.where((dofs_f_V[:,1] == 0.19) & (dofs_f_V[:,0] >= 0.24) & (dofs_f_V[:,0] <= 0.6))[0]
		i_m_V_fsi2 = np.where((self.dofs_m_V[:,1] == 0.19) & (self.dofs_m_V[:,0] >= 0.24) & (self.dofs_m_V[:,0] <= 0.6))[0]

		# end
		i_s_V_fsi3 = np.where((self.dofs_s_V[:,0] == 0.6) & (self.dofs_s_V[:,1] >= 0.19 +DOLFIN_EPS) & (self.dofs_s_V[:,1] <= 0.21 - DOLFIN_EPS))[0]
		i_f_V_fsi3 = np.where((dofs_f_V[:,0] == 0.6) & (dofs_f_V[:,1] >= 0.19+DOLFIN_EPS) & (dofs_f_V[:,1] <= 0.21 -DOLFIN_EPS))[0]
		i_m_V_fsi3 = np.where((self.dofs_m_V[:,0] == 0.6) & (self.dofs_m_V[:,1] >= 0.19+DOLFIN_EPS) & (self.dofs_m_V[:,1] <= 0.21 -DOLFIN_EPS))[0]

		i_f_V_in = np.where((dofs_f_V[:,0] == 0.0))[0]


		# Concatenate
		self.i_s_V_fsi = np.concatenate((i_s_V_fsi2,i_s_V_fsi3,i_s_V_fsi1),axis = 0)
		i_f_V_fsi = np.concatenate((i_f_V_fsi2,i_f_V_fsi3,i_f_V_fsi1),axis = 0)
		i_m_V_fsi = np.concatenate((i_m_V_fsi2,i_m_V_fsi3,i_m_V_fsi1),axis = 0)

		# This is a complete set of indices. Find ones that match between fluid and structure:

		# This might not be the best approach. For instance, interpolating the structure solution.

		# Fluid has more nodes because for lagrange elements of order greater than 2 there are
		# nodes that do not correspond to vertices.

		# overlapping nodes between structure and fluid. Vertices match, but nodes do not because higher order.
		self.i_f_V_fsi_com = []
		for i_s in range(self.dofs_s_V[self.i_s_V_fsi].shape[0]/2):
			for i_f in range(dofs_f_V[i_f_V_fsi].shape[0]/2):
				if self.dofs_s_V[self.i_s_V_fsi[2*i_s]][0] == dofs_f_V[i_f_V_fsi[2*i_f]][0] and self.dofs_s_V[self.i_s_V_fsi[2*i_s]][1] == dofs_f_V[i_f_V_fsi[2*i_f]][1]:
					self.i_f_V_fsi_com.append(i_f_V_fsi[2*i_f])
					self.i_f_V_fsi_com.append(i_f_V_fsi[2*i_f+1])

		self.i_m_V_fsi_com = []
		for i_s in range(self.dofs_s_V[self.i_s_V_fsi].shape[0]/2):
			for i_m in range(self.dofs_m_V[i_m_V_fsi].shape[0]/2):
				if self.dofs_s_V[self.i_s_V_fsi[2*i_s]][0] == self.dofs_m_V[i_m_V_fsi[2*i_m]][0] and self.dofs_s_V[self.i_s_V_fsi[2*i_s]][1] == self.dofs_m_V[i_m_V_fsi[2*i_m]][1]:
					self.i_m_V_fsi_com.append(i_m_V_fsi[2*i_m])
					self.i_m_V_fsi_com.append(i_m_V_fsi[2*i_m+1])



	def saveResults(self, structureSolver, fluidSolver, meshSolver):
		# Save fluid velocity and pressure
		fluidSolver.u_res.assign(fluidSolver.u1)
		self.file_v_f << (fluidSolver.u_res, self.t)
		fluidSolver.p_res.assign(fluidSolver.p1)
		self.file_p_f << (fluidSolver.p_res, self.t)

		meshSolver.u_res.assign(meshSolver.meshVelocity)
		self.file_v_m << (meshSolver.u_res, self.t)

		#Save structure displacement and velocity
		# extract displacements and velocities for results...
		u2, v2 = structureSolver.U.split()

		# maybe like this, or maybe save U.
		self.file_u_s << (u2, self.t)
		self.file_v_s << (v2, self.t)
		#structureSolver.u_res.assign(structureSolver.u)
		#self.file_u_s << (structureSolver.u_res, self.t)
		#structureSolver.v_res.assign(structureSolver.v)
		#self.file_v_s << (structureSolver.v_res, self.t)
