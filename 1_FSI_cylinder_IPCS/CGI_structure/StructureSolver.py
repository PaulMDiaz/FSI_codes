# Structure solver. Intended to eventually be generic, but for now addresses St venant model only.
#
#
#  Inputs:
#    IOInfo - An object which contains all inputs and outputs for the
#    Linear Elastic eqution solver implemented herein

from dolfin import *
from kinematics import *
import numpy as np


class Structure_Solver:
	"""Object for FEM structure solver"""
	def __init__(self, Mesh, ElementType, ElementDegree, StructureSolverMethod, StructureBodyForce):

		# Initialize mesh and FEM solver parameter
		self.mesh = Mesh # static mesh
		self.ElementType = ElementType
		self.ElementDegree = ElementDegree
		self.solverMethod = StructureSolverMethod
		self.B = StructureBodyForce

		# Variational spaces
		self.S_space = FunctionSpace(self.mesh, self.ElementType, self.ElementDegree)
		self.V_space = VectorFunctionSpace(self.mesh, self.ElementType, self.ElementDegree)
		self.T_space = TensorFunctionSpace(self.mesh, self.ElementType, self.ElementDegree)

		self.V_element = VectorElement(self.ElementType, self.mesh.ufl_cell() , self.ElementDegree)

		# Mixed space for computing both displacement and rate of change of displacement
		self.V_element = VectorElement(self.ElementType, self.mesh.ufl_cell() , self.ElementDegree)
		#self.V2_temp = self.V_element*self.V_element

		#self.mixed_space = FunctionSpace(self.mesh, self.V2_temp)
		self.mixed_space = FunctionSpace(self.mesh, MixedElement([self.V_element, self.V_element]))
		self.mixed_space = FunctionSpace(self.mesh, self.V_element*self.V_element)

		self.V = TestFunction(self.mixed_space)
		self.dU = TrialFunction(self.mixed_space)
		self.U = Function(self.mixed_space)
		self.U0 = Function(self.mixed_space)

		# Load initial conditions to u0 and v0. Otherwise set to 0.
		self.u0 = Constant((0,)*self.V_space.mesh().geometry().dim())
		self.v0 = Constant((0,)*self.V_space.mesh().geometry().dim())

		# Functions for solver
		self.u_t, self.v_t = split(self.V) 	# Test functions
		self.u, self.v = split(self.U)		# Functions

		self.cells = CellFunction("size_t", self.mesh)
		self.dx = Measure('dx', domain = self.mesh, subdomain_data = self.cells)

		# Project u0 and v0 into U0
		self.a_proj = inner(self.dU, self.V)*self.dx
		self.L_proj = inner(self.u0, self.u_t)*self.dx + inner(self.v0, self.v_t)*self.dx
		solve(self.a_proj == self.L_proj, self.U0)

		# load data
		#U_temp = np.loadtxt('nodal_U_64')
		#self.U0.vector()[:] = U_temp

		self.u0, self.v0 = split(self.U0)

		## Function for results
		self.u_res = Function(self.V_space, name = 'u')
		self.v_res = Function(self.V_space, name = 'v')

	def _construct_local_kinematics(self, u):
		self.I = SecondOrderIdentity(u)
		self.epsilon = InfinitesimalStrain(u)
		self.F = DeformationGradient(u)
		self.J = Jacobian(u)
		self.C = RightCauchyGreen(u)
		self.E = GreenLagrangeStrain(u)
		self.b = LeftCauchyGreen(u)
		self.e = EulerAlmansiStrain(u)
		[self.I1, self.I2, self.I3] = CauchyGreenInvariants(u)
		[self.I1bar, self.I2bar] = IsochoricCauchyGreenInvariants(u)

	# appropriate for St venant only:
	def SecondPiolaKirchhoffStress(self, u, p_s):
		self._construct_local_kinematics(u)
		psi = p_s.lambda_s/2*(tr(self.E)**2)+p_s.mu_s*tr(self.E*self.E)
		E = self.E
		S = diff(psi, E)
		return S

	def FirstPiolaKirchhoffStress(self, u, p_s):
		S = self.SecondPiolaKirchhoffStress(u, p_s)
		F = self.F
		P = F*S
		return P

###########################
	#	self.d00_s = Function(self.V_space)
###########################
	def Structure_Problem_Solver(self, p_s, F):

		if self.solverMethod == "Linear":
			self.Linear_Elastic_Solver(p_s, F)
		elif self.solverMethod == "NeoHookean":
			#self.Incompressible_NeoHookean_Solver(p_s, F)
			self.Compressible_NeoHookean_Solver(p_s, F)
		elif self.solverMethod =="St_venant":
			self.Compressible_St_venant(p_s, F)

		else:
			print "Error. The only solvers available for the structure are Linear or NeoHookean"

	def Compressible_St_venant(self, p_s, F):

		print ""
		print ""
		print "ENTERING STRUCTURE St Venant SOLVER"
		print ''
		print ''

		# displacements and velocities at mid points
		self.u_mid = 0.5*(self.u0 + self.u)
		self.v_mid = 0.5*(self.v0 + self.v)

		#I = Identity(p_s.Dim)

		# Project stress from fluid tensor space onto structure tensor space
		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )
		#self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			#form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		# Traction on boundary

		# solves and gives 0 displacment
		#self.T_hat = Constant((0.0, 0.0))

		# Solves and gives too great a displacement.
		n_s = FacetNormal(self.mesh)

		self.u1, self.v1 = self.U.split()

		# Calculate necessary kinematics
		# Jacobian
		self.J = Jacobian(self.u1)
		# Deformation gradient tensor
		self.F = DeformationGradient(self.u1)

		self.T_hat = (self.J*inv(self.F)*self.sigma_FSI).T * n_s
		# gives same results as index method
		#self.T_hat = Constant((0.0, 0.0))
		#self.T_hat = as_tensor( self.J*inv(self.F)[k, j]*self.sigma_FSI[j, i]*n_s[k] , (i, ) )

		#self.Dim = self.mesh.topology().dim()
		#i, j, k, l, m = indices(5)
		#self.delta = Identity(self.Dim)

		# The variational form corresponding to hyperelasticity

		self.k = Constant(p_s.dt)
		self.P1 = self.FirstPiolaKirchhoffStress(self.u_mid, p_s)

		self.L = p_s.rho_s*inner(self.v - self.v0, self.u_t)*self.dx \
		+ self.k*inner(self.P1, grad(self.u_t))*self.dx \
		- self.k*inner(self.B, self.u_t)*self.dx \
		+ inner(self.u - self.u0, self.v_t)*self.dx \
		- self.k*inner(self.v_mid, self.v_t)*self.dx

		# Implement Neuman BCs
		self.L = self.L-self.k*inner(self.T_hat, self.u_t)*self.ds(2)

		# introduced in line 499 of solution_algorithms
		self.a = derivative(self.L, self.U, self.dU)

		# Setup problem

		problem = NonlinearVariationalProblem(self.L, self.U, self.bcs, self.a)
		solver = NonlinearVariationalSolver(problem)
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["maximum_iterations"] = 100

		solver.solve()

		#self.u, self.v = self.U.split()

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"
