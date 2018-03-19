# Structure solver. Intended to eventually be generic, but for now addresses St venant model only.
#
#
#  Inputs:
#    IOInfo - An object which contains all inputs and outputs for the
#    Linear Elastic eqution solver implemented herein

from dolfin import *
from kinematics import *
import numpy as np


class StructureSolver:
	"""Object for FEM structure solver"""
	def __init__(self, elementType, elementDegree, structureSolverMethod, structureBodyForce, p_s):

		# Initialize mesh and FEM solver parameter
		self.mesh = p_s.meshStructure # static mesh
		self.elementType = elementType
		self.elementDegree = elementDegree
		self.solverMethod = structureSolverMethod

		self.B = structureBodyForce

		# Variational spaces
		self.scalarSpace = FunctionSpace(self.mesh, self.elementType, self.elementDegree)
		self.vectorSpace = VectorFunctionSpace(self.mesh, self.elementType, self.elementDegree)
		self.tensorSpace = TensorFunctionSpace(self.mesh, self.elementType, self.elementDegree)

		# Mixed space for computing both displacement and rate of change of displacement
		self.vectorElement = VectorElement(self.elementType, self.mesh.ufl_cell() , 1)
		self.mixedSpace = FunctionSpace(self.mesh, MixedElement([self.vectorElement, self.vectorElement]))

		self.V = TestFunction(self.mixedSpace)
		self.dU = TrialFunction(self.mixedSpace)
		self.U = Function(self.mixedSpace)
		self.U0 = Function(self.mixedSpace)

		# Load initial conditions to u0 and v0. Otherwise set to 0.
		self.u0 = Constant((0,)*self.vectorSpace.mesh().geometry().dim())
		self.v0 = Constant((0,)*self.vectorSpace.mesh().geometry().dim())

		# Functions for solver
		self.xi, self.eta = split(self.V) 	# Test functions
		self.u, self.v = split(self.U)		# Functions

		self.cells = CellFunction("size_t", self.mesh)
		self.dx = Measure('dx', domain = self.mesh, subdomain_data = self.cells)

		# Project u0 and v0 into U0
		self.a_proj = inner(self.dU, self.V)*self.dx
		self.L_proj = inner(self.u0, self.xi)*self.dx + inner(self.v0, self.eta)*self.dx
		solve(self.a_proj == self.L_proj, self.U0)

		# load saved data
		# U_temp = np.loadtxt('Restart_FSI/nodal_Ustr_1') 
		# self.U0.vector()[:] = U_temp

		self.u0, self.v0 = split(self.U0)

		# vector space for traction, just use vectorSpace
		self.traction = Function(self.vectorSpace)

		self.u_vel = Function(self.vectorSpace)

		# self.interfaceStructureTensorSpace = TensorFunctionSpace(p_s.interfaceStructure, "CG", 1)
		# interfaceStructureTensorCoords = self.interfaceStructureTensorSpace.tabulate_dof_coordinates().reshape((self.interfaceStructureTensorSpace.dim(),-1))

		# these coordinates match the fluid tensor space ones.

		# Function for fluid interface velocity
		# self.structureInterfaceStress = Function(self.interfaceStructureTensorSpace)
		# self.structureInterfaceStress.vector().zero()

# structureSolver.interfaceStructureTensorSpace = TensorFunctionSpace(p_s.interfaceStructure, "CG", 1)
# interfaceStructureTensorCoords = structureSolver.interfaceStructureTensorSpace.tabulate_dof_coordinates().reshape((structureSolver.interfaceStructureTensorSpace.dim(),-1))
#
# # Function for fluid interface velocity
# structureSolver.structureInterfaceStress = Function(structureSolver.interfaceStructureTensorSpace)
# structureSolver.structureInterfaceStress.vector().zero()

		## Function for results
		self.u_res = Function(self.vectorSpace, name = 'u')
		self.v_res = Function(self.vectorSpace, name = 'v')

	def _construct_local_kinematics(self, u):
		self.I = secondOrderIdentity(u)
		self.epsilon = infinitesimalStrain(u)
		self.F = deformationGradient(u)
		self.J = jacobian(u)
		self.C = rightCauchyGreen(u)
		self.E = greenLagrangeStrain(u)
		self.b = leftCauchyGreen(u)
		self.e = EulerAlmansiStrain(u)

		[self.I1, self.I2, self.I3] = cauchyGreenInvariants(u)
		[self.I1bar, self.I2bar] = isochoricCauchyGreenInvariants(u)

	# appropriate for St venant only:
	def secondPiolaKirchhoffStress(self, u, p_s):
		self._construct_local_kinematics(u)
		psi = p_s.lambda_s/2*(tr(self.E)**2)+p_s.mu_s*tr(self.E*self.E)
		E = self.E
		S = diff(psi, E)
		return

	def firstPiolaKirchhoffStress(self, u, p_s):
		S = self.secondPiolaKirchhoffStress(u, p_s)
		F = self.F
		P = F*S
		return P

###########################
	#	self.d00_s = Function(self.vectorSpace)
###########################
	def structureProblemSolver(self, p_s, fluidSolver):

		if self.solverMethod == "Linear":
			self.linearElasticSolver(p_s, fluidSolver)
		elif self.solverMethod == "NeoHookean":
			#self.Incompressible_NeoHookean_Solver(p_s, fluidSolver)
			self.compressibleNeoHookeanSolver(p_s, fluidSolver)
		elif self.solverMethod =="St_venant":
			self.compressibleStVenant(p_s, fluidSolver)

		else:
			print "Error. The only solvers available for the structure are Linear or NeoHookean"

	def compressibleStVenant(self, p_s, fluidSolver):

		print ""
		print ""
		print "ENTERING STRUCTURE St Venant SOLVER"
		print ''
		print ''

		self.u0, self.v0 = split(self.U0)
        #
		# # Project stress from reference fluid tensor space onto structure tensor space
		# not confident that this will yield correct forces... compute a negative perhaps?

		# in cbc, it appears that mesh is solved for as displacement...
		#self.sigma_FSI = project(fluidSolver.sigma_FSI, self.tensorSpace, solver_type = "mumps",\
			#form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )

		# # copy stress values on interface from fluid to structure.
		# self.structureInterfaceStress.vector()[:] = fluidSolver.refInterfaceStress.vector().get_local()

# copy values from fluid tensor space to structure tensor space... little difficult...

		# Kinematics
		# displacements and velocities at mid points
		self.u_mid = 0.5*(self.u0 + self.u)
		self.v_mid = 0.5*(self.v0 + self.v)


		self.d = self.u.geometric_dimension()
		self.I = Identity(self.d)         #identity
		self.F = self.I + grad(self.u_mid)     #deformation gradient
		self.C = self.F.T*self.F               #Right Cauchy-Green tensor
		E = (self.C - self.I)/2           #Green-Lagrange strain tensor
		self.E = variable(E)

		# St venant model for stored strain energy
		self.psi = p_s.lambda_s/2*(tr(self.E)**2) + p_s.mu_s*tr(self.E*self.E)
		# Second Piola-Kirchhoff stress
		self.S = diff(self.psi, self.E)
		# First Piola-Kirchhoff stress
		self.P1 = self.F*self.S
		self.k = Constant(p_s.dt)
		# Solves and gives too great a displacement.
		self.n_s = FacetNormal(self.mesh)

		# traction on structure.
		# believe that traction has already been negated in fluid problem.
		# therefore no need to do it here.
		class FluidInterfaceTractionExpression(Expression):
		    def eval(self, values, x):
		        try:
		            values[:] = fluidSolver.traction(x)
		        except:
		            values[:] = 0
		    def value_shape(self):
		        return (2,)

		self.fluidInterfaceTractionExpression = FluidInterfaceTractionExpression(degree=2)

#		# Compute integral of transferred projection
#		form = dot(self.G_S, self.N_S)*d_FSI
#		integral_2 = assemble(form, exterior_facet_domains=self.problem.fsi_boundary_S)

		self.L = p_s.rho_s*inner(self.v - self.v0, self.xi)*self.dx \
		+ self.k*inner(self.P1, grad(self.xi))*self.dx \
		- self.k*inner(self.B, self.xi)*self.dx \
		+ inner(self.u - self.u0, self.eta)*self.dx \
		- self.k*inner(self.v_mid, self.eta)*self.dx

		# Implement Neuman BCs
		self.L = self.L - self.k*inner(self.fluidInterfaceTractionExpression, self.xi)*self.ds(20)

		# introduced in line 499 of solution_algorithms
		self.a = derivative(self.L, self.U, self.dU)

		# Setup problem
		problem = NonlinearVariationalProblem(self.L, self.U, self.bcs, self.a)
		solver = NonlinearVariationalSolver(problem)
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-8
		solver.parameters["newton_solver"]["relative_tolerance"] = 1e-8
		solver.parameters["newton_solver"]["maximum_iterations"] = 1000

		solver.solve()



		# Extract velocity to compare
		self.u12, self.v12 = self.U.split(deepcopy = True)
		self.u02, self.v02 = self.U0.split(deepcopy = True)

# 		self.u12.set_allow_extrapolation(deepcopy = True)
# structureSolver.u12.set_allow_extrapolation(deepcopy = True)
# # structureSolver.u12
# u11, v11 = structureSolver.U.split(deepcopy = True)
		self.u_vel.vector()[:] = (1.0/self.k)*(self.u12.vector().get_local() - self.u02.vector().get_local())

		#self.u, self.v = self.U.split()

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"
