#  Configures variational problem and boundary conditions for the linear elastic structural displacement equations
#  to solve for the structural displacement and forces and mesh displacement in
#  a Driven Cavity Flow
#
#  Created: 12 May 2015
#  Modified: 3 June 2015
#
#  Inputs:
#    IOInfo - An object which contains all inputs and outputs for the
#    Linear Elastic eqution solver implemented herein

from dolfin import *
from kinematics import *
import numpy as np


class Structure_Solver:
	"""Object for FEM structure solver"""
	def __init__(self, Mesh, ElementType, ElementDegree, StructureSolverMethod, StructureBodyForce, p_s, F):

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

		# Load initial conditions to u0 and v0. Otherwise set to 0.
		self.u0 = Constant((0,)*self.V_space.mesh().geometry().dim())
		self.v0 = Constant((0,)*self.V_space.mesh().geometry().dim())

		# Functions for solver
		self.v = TestFunction(self.V_space)
		self.u1 = Function(self.V_space)
		self.v1 = Function(self.V_space)
		self.a1 = Function(self.V_space)
		self.du = TrialFunction(self.V_space)

		# Initial displacement and velocity
		self.u0 = interpolate(self.u0, self.V_space)
		self.v0 = interpolate(self.v0, self.V_space)
		self.v1 = interpolate(self.v0, self.V_space)

		# Parameters for HHT time integration
		# the current settings ensure 2nd order accurate, stable for linear problems and introduces no numerical dissipation
		#self.alpha = 1.0
		self.beta = 0.25
		self.gamma = 0.5

		# Deterimine initial acceleration
		self.a0 = TrialFunction(self.V_space)

		self.cells = CellFunction("size_t", self.mesh)
		self.dx = Measure('dx', domain = self.mesh, subdomain_data = self.cells)

		## Function for results
		self.u_res = Function(self.V_space, name = 'u')
		self.v_res = Function(self.V_space, name = 'v')

		self.P0 = self.FirstPiolaKirchhoffStress(self.u0, p_s)
		# Deterimine initial acceleration
		self.a_accn = inner(self.a0, self.v)*self.dx
		self.L_accn = - inner(self.P0, grad(self.v))*self.dx + inner(self.B, self.v)*self.dx

		# Solves and gives too great a displacement.


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
	def Structure_Problem_Solver(self, p_s, F):
		if self.solverMethod =="St_venant":
			self.Compressible_St_venant(p_s, F)
		else:
			print "Error. The only solvers available for the structure are Linear or NeoHookean"

	def Compressible_St_venant(self, p_s, F):

		print ""
		print ""
		print "ENTERING STRUCTURE St Venant SOLVER"
		print ''
		print ''

		# Project stress from fluid tensor space onto structure tensor space
		#self.T_hat = Constant((0.0, 0.0))

		n_s = FacetNormal(self.mesh)

		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		self._construct_local_kinematics(self.u1)
		# To be consistent, use u_mid for T_hat too...
		self.T_hat = (self.J*inv(self.F)*self.sigma_FSI).T * n_s

	#	self.T_hat = Constant((0.0, 0.0))

		self.L_accn = self.L_accn+inner(self.T_hat, self.v)*self.ds(2)

		self.k = Constant(p_s.dt)
		self.a1 = self.a0*(1.0-1.0/(2*self.beta)) - (self.u0-self.u1 +self.k*self.v0)/(self.beta*self.k**2)

		self.P = self.FirstPiolaKirchhoffStress(self.u1, p_s)

#         #  A general version of the trick below is what should
#         # be used instead. The commentend-out lines only work well for
#         # quadratically nonlinear models, e.g. St. Venant Kirchhoff.
	##
 		#self.S0 = self.SecondPiolaKirchhoffStress(self.u0, p_s)
 		#self.S1 = self.SecondPiolaKirchhoffStress(self.u1, p_s)
 		#self.Sm = 0.5*(self.S0 + self.S1)
		#self._construct_local_kinematics(0.5*(self.u0 + self.u1))
		#self.P  = self.F*self.Sm
	##
		# The variational form corresponding to hyperelasticity (did have int(problem.is_dynamic()), ie 1 when dynamic multiplying first term)
		self.L = p_s.rho_s*inner(self.a1, self.v)*self.dx \
			+ inner(self.P, grad(self.v))*self.dx - inner(self.B, self.v)*self.dx

		self.L = self.L-inner(self.T_hat, self.v)*self.ds(2)

		self.a = derivative(self.L, self.u1, self.du)

		problem = NonlinearVariationalProblem(self.L, self.u1, self.bcs, self.a)
		solver = NonlinearVariationalSolver(problem)
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["maximum_iterations"] = 100

		solver.solve()

		# Update acceleration and velocity. Further updates completed in FSI_Driver.
		self.a1 = self.a0*(1.0 - 1.0/(2*self.beta)) \
		    - (self.u0 - self.u1 + self.k*self.v0)/(self.beta*self.k**2)
		self.a1 = project(a1, self.vector)
		self.v1 = self.v0 + self.k*((1 - self.gamma)*self.a1 + self.gamma*self.a0)
		self.v1 = project(v1, self.vector)

		# First directional derivative of Pi about d in the direction of v
		#Form_s = derivative(self.Pi, self.d_, self.v)

		# Jacobian of the directional derivative Fd
		#Gain_s = derivative(Form_s, self.d_, self.du)

		#begin("Computing structure displacement")

		#solve(Form_s == 0, self.d_, self.bcs, J = Gain_s, \
			#solver_parameters  = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
			#form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
		#end()

		#print ""
		#print ""
		#print "EXITING STRUCTURE SOLVER"

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"
