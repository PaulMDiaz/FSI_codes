#  This class contains the functions necessary to setup the FSI problem
#  and to handle the transfer of information between the fluid and structure
#  solvers during the time-stepping iteration.

#  Created: 28 March 2014
#  Modified: 23 August 2015
#from dolfin import *
from fenics import *
import numpy as np
import pylab as plt

class MeshSolver(object):

	def __init__(self, meshElementType, meshElementDegree, p_s):

		self.mesh = p_s.meshRef	# reference mesh
		self.elementType = meshElementType
		self.elementDegree = meshElementDegree


		self.tensorSpace = TensorFunctionSpace(self.mesh, self.elementType, self.elementDegree)
		self.vectorSpace = VectorFunctionSpace(self.mesh, self.elementType, self.elementDegree)
		self.scalarSpace = FunctionSpace(self.mesh, self.elementType, self.elementDegree)

		# time-dependent linear elasticity

		# move to FSI_Driver eventually.
		self.mu = 1.0
		self.lmbda = 1.0
		self.alpha = 1.0

		self.v = TestFunction(self.vectorSpace)
		self.u = TrialFunction(self.vectorSpace)
		self.u0 = Function(self.vectorSpace)	# previous time step mesh displacement
		self.u1 = Function(self.vectorSpace)	# current time step mesh displacement

		# u_temp = np.loadtxt('Restart_FSI/nodal_m_disp1')
		# self.u1.vector()[:] = u_temp
		# u_temp = np.loadtxt('Restart_FSI/nodal_m_disp0')
		# self.u0.vector()[:] = u_temp

		self.meshVelocity = Function(self.vectorSpace)
		self.meshDisplacment = Function(self.vectorSpace)

		self.num_dofs = self.vectorSpace.dim()

		self.u_res = Function(self.vectorSpace, name = 'u')
		# self.v_res = Function(self.vectorSpace, name = 'v')
		# Define boundary conditions: addressed in Cylinder.
		# Will have to measure displacement later.

		def sigma_M(U_M, mu_M, lmbda_M):
			"Return mesh stress in reference domain"
			I = variable(Identity(U_M.geometric_dimension()))
			return 2*mu_M*sym(grad(U_M)) + lmbda_M*tr(grad(U_M))*I

		self.cells = CellFunction("size_t", self.mesh)
		self.facets = FacetFunction("size_t", self.mesh)

		self.dx = Measure('dx', domain = self.mesh, subdomain_data = self.cells)
		self.ds = Measure('ds', domain = self.mesh, subdomain_data = self.facets)

		# Define cG(1) scheme for time-stepping
		self.k = Constant(0)
		self.a = self.alpha*inner(self.v, self.u)*self.dx \
		 + 0.5*self.k*inner(sym(grad(self.v)), sigma_M(self.u, self.mu, self.lmbda))*self.dx
		self.L = self.alpha*inner(self.v, self.u0)*self.dx \
		 - 0.5*self.k*inner(sym(grad(self.v)), sigma_M(self.u0, self.mu, self.lmbda))*self.dx

		# Add right-hand side: for body forces. No need just now.
		# F_M = self.problem.mesh_right_hand_side()
		# L += k*inner(v, F_M)*dx

		# self.meshVelocity = Function(self.vectorSpace) # Current time step
		# self.meshDisplacement1 = Function(self.vectorSpace) # Current time step
		# self.meshDisplacement1.vector().zero()
		# self.meshDisplacement0 = Function(self.vectorSpace) # previos time step
		# self.meshDisplacement0.vector().zero()
		# self.meshDisplacement = Function(self.vectorSpace) # previos time step
		# self.meshDisplacement.vector().zero()
        #
		# # self.v_mesh = Function(self.vectorSpace) # for projection of structure velocity
        #
		# # Trial and test functions for mesh motion
		# self.del_u_mesh = TestFunction(self.vectorSpace)
		# self.du_mesh = TrialFunction(self.vectorSpace)

		# self.interfaceMeshVectorSpace = VectorFunctionSpace(p_s.interfaceMesh, "CG", 1)
		# self.meshInterfaceVelocity = interpolate(self.meshVelocity,self.interfaceMeshVectorSpace)
		# self.interfaceMeshCoords = self.interfaceMeshVectorSpace.tabulate_dof_coordinates().reshape((self.interfaceMeshVectorSpace.dim(),-1))

	def meshProblemSolver(self, structureSolver, fluidSolver, p_s):
		# def Move_Mesh(self, Structure, Fluid):
		# Description: Moves structure and fluid meshes given structure displacement
		# Inputs: Structure - structure solver object
		#		  Fluid - fluid solver object

		print ''
		print ''
		print 'ENTERING MESH MOVEMENT'
		print ''
		print ''

		self.k.assign(p_s.dt)

		# Assemble linear system and apply boundary conditions
		self.A = assemble(self.a)
		self.b = assemble(self.L)

		# this is how it reads in cbc. In my other solver I have bc in solve.
		# self.bcMesh.apply(self.A, self.b)

		[bc.apply(self.A) for bc in self.bcMesh]
		[bc.apply(self.b) for bc in self.bcMesh]

		# Compute solution
		solve(self.A, self.u1.vector(), self.b)

		# u1 is mesh displacement
		# require a mesh velocity to apply to fluid.

		# Compute mesh velocity
		self.meshVelocity.vector()[:] = (1.0/self.k)*(self.u1.vector().get_local()-self.u0.vector().get_local())
		self.meshVelocity.set_allow_extrapolation(True)


		# # artificial viscosity
		# a = 1.0e-11		# MPa/s
        #
		# self.Form_m = a*sym(grad(self.meshVelocity))[i, j]*self.del_u_mesh[i].dx(j)*self.dx
		# self.Gain_m = derivative(self.Form_m, self.meshVelocity, self.du_mesh)
        #
		# solve(self.Form_m == 0, self.meshVelocity, self.bcMesh, J = self.Gain_m,
		# 	solver_parameters = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-10} }, \
		# 	form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )
        #
		# # try with no mesh movement to trouble shoot fluid solver.
        #
		# interpolate mesh velocity solution to FSI interface
		# self.meshInterfaceVelocity = interpolate(self.meshVelocity,self.interfaceMeshVectorSpace)

		# do I need these? Doubtful...
		# interpolate mesh displacement to FSI interface to update interface mesh.
		# self.meshInterfaceDisplacement = interpolate(self.u1,self.interfaceMeshVectorSpace)

		self.meshDisplacment = 0.5*(self.u1+self.u0)

	def Generate_Files(Structure, Fluid):
		File("fluid_pressure.pvd") << Fluid.p1
		File("fluid_velocity.pvd") << Fluid.u1
		File("fluid_mesh.pvd") << Fluid.mesh
		File("mesh.pvd") << mesh
		File("structure_mesh.pvd") << Structure.mesh
		File("structure_displacement.pvd") << Structure.d
		print ""
		print "Simulation Finished"
		print "Close plots to exit"
