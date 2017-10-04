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
import numpy as np


class Structure_Solver:
	"""Object for FEM structure solver"""
	def __init__(self, Mesh, ElementType, ElementDegree, StructureSolver, StructureBodyForce):

		# Initialize mesh and FEM solver parameter
		self.mesh = Mesh # static mesh
		self.ElementType = ElementType
		self.ElementDegree = ElementDegree
		self.solver = StructureSolver
		self.B = StructureBodyForce

		# Variational spaces
		self.S_space = FunctionSpace(self.mesh, self.ElementType, self.ElementDegree)
		self.V_space = VectorFunctionSpace(self.mesh, self.ElementType, self.ElementDegree)
		self.T_space = TensorFunctionSpace(self.mesh, self.ElementType, self.ElementDegree)

		# Function for results
		self.d_ = Function(self.V_space, name = 'd')

		# Functions for solver
		self.d = Function(self.V_space)
		self.d0 = Function(self.V_space)
		self.u = TrialFunction(self.V_space)
		self.du = TestFunction(self.V_space)

	def Structure_Problem_Solver(self, DC, F):

		if self.solver == "Linear":
			self.Linear_Elastic_Solver(DC, F)
		elif self.solver == "NeoHookean":
			#self.Incompressible_NeoHookean_Solver(DC, F)
			self.Compressible_NeoHookean_Solver(DC, F)
		else:
			print "Error. The only solvers available for the structure are Linear or NeoHookean"

		self.d_dot = project( (self.d - self.d0)/DC.dt, self.V_space, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

	def Linear_Elastic_Solver(self, DC, F):

		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		I = Identity(DC.Dim)

		self.sigma = 2.0*DC.mu_s*sym(grad(self.u)) + DC.lambda_s*tr(sym(grad(self.u)))*I

		self.F = I + grad(self.u)
		self.J = det(self.F)
		self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N
		#self.T_hat = Constant((0.0,0.0))

		print ""
		print ""
		print "ENTERING STRUCTURE LINEAR SOLVER"
		print ''
		print ''

		# Setup variational problem
		self.a = inner(self.sigma, grad(self.du))*self.dV
		self.L = inner(self.T_hat, self.du)*self.dA(3) # ds(3) = FSI boundary

		# Solve variational problem
		begin("Computing structure displacement")
		solve(self.a == self.L, self.d, self.bcs, solver_parameters={"symmetric":True})
		end()
			#self.Incompressible_NeoHookean_Solver(DC, F)

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"

	def Compressible_NeoHookean_Solver(self, DC, F):

		print ""
		print ""
		print "ENTERING STRUCTURE NEOHOOKEAN SOLVER"
		print ''
		print ''

		I = Identity(DC.Dim)

		# Project stress from fluid tensor space onto structure tensor space
		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )


		# Kinematics
		self.F = I + grad(self.d)			# Deformation gradient
		self.C = self.F.T*self.F			# Right Cauchy-Green tensor

		# Invariants of deformation tensor
		self.Ic = tr(self.C)
		self.J = det(self.F)

		# Elasticity parameters E and nu should be defined

		# Stored strain energy density (compressible neo-Hookean model)
		self.psi = (DC.mu_s/2)*(self.Ic - 3) - DC.mu_s*ln(self.J) + (DC.lambda_s/2)*(ln(self.J))**2

		#original line
		#self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N

		# Traction force (per unit reference area)
		# (first piola kirchoff stress tensor, also called called lagrangian stress tensor)
		#self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N
		self.T_fat = self.J*inv(self.F)*self.sigma_FSI*self.N

		# print out pressure/ T_hat on interface:
		# Structure array of coordinates
		dofs_s_S = self.S_space.tabulate_dof_coordinates().reshape((self.S_space.dim(),-1))
		dofs_s_V = self.V_space.tabulate_dof_coordinates().reshape((self.V_space.dim(),-1))
		dofs_s_T = self.T_space.tabulate_dof_coordinates().reshape((self.T_space.dim(),-1))

		i_s_S = np.where((dofs_s_S[:,1] == DC.h))[0] #  & (x <= 0.5)
		i_s_V = np.where((dofs_s_V[:,1] == DC.h))[0] #  & (x <= 0.5)
		i_s_T = np.where((dofs_s_T[:,1] == DC.h))[0] #  & (x <= 0.5)

		# Extract dofdofs_f_S indices for fluid velocity on top plate
		i_s_V_FSI = np.where((dofs_s_V[:,1] == DC.h))[0]
		i_s_S_FSI = np.where((dofs_s_S[:,1] == DC.h))[0]
		i_s_T_FSI = np.where((dofs_s_T[:,1] == DC.h))[0]


		#self.pressure_s = self.T_hat.vector()[i_s_S_FSI]

		#print 'Structure T_hat ', self.T_hat

		#self.p_test = dot(self.T_hat, self.d)*self.dA(3)

		#print 'pressure_FS =', self.p_test.vector()

		self.T_hat = Constant((0.0,0.0))

		# Total potential energy dA(3) integrate on interface only.
		self.Pi = self.psi*self.dV - dot(self.T_hat, self.d)*self.dA(3) - dot(self.B, self.d)*self.dV

		# First directional derivative of Pi about d in the direction of v
		Form_s = derivative(self.Pi, self.d, self.du)

		# Jacobian of the directional derivative Fd
		Gain_s = derivative(Form_s, self.d, self.u)

		begin("Computing structure displacement")
		# Solve variational problem. Why these solver parameters?

		#Could use:
		parameters["form_compiler"]["cpp_optimize"] = True
		ffc_options = {"optimize": True}

		#solve(Form_s == 0, self.d, self.bcs, J=Gain_s,
		#      form_compiler_parameters=ffc_options)

		solve(Form_s == 0, self.d, self.bcs, J = Gain_s, \
			solver_parameters  = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
		end()

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"

	def Incompressible_NeoHookean_Solver(self, DC, F):

		print ""
		print ""
		print "ENTERING STRUCTURE NEOHOOKEAN SOLVER"
		print ''
		print ''

		I = Identity(DC.Dim)

		# Project stress from fluid tensor space onto structure tensor space.
		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		# Kinematics
		self.F = I + grad(self.d)			# Deformation gradient
		self.C = self.F.T*self.F			# Right Cauchy-Green tensor

		# Invariants of deformation tensor
		self.Ic = tr(self.C)
		self.J = det(self.F)

		# Elasticity parameters E and nu should be defined.

		# Stored strain energy density
		self.psi = (DC.mu_s/2)*(self.Ic - 3) - DC.mu_s*ln(self.J) + (DC.lambda_s/2)*(ln(self.J))**2
		#self.psi = (DC.mu_s/2)*(self.Ic-3)

		#self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N
		self.T_hat = Constant((0.0,0.0))
		# Total potential energy
		self.Pi = self.psi*self.dV - dot(self.T_hat, self.d)*self.dA(3) - dot(self.B, self.d)*self.dV

		# First directional derivative of Pi about d in the direction of v
		Form_s = derivative(self.Pi, self.d, self.du)

		# Jacobian of the directional derivative Fd
		Gain_s = derivative(Form_s, self.d, self.u)

		begin("Computing structure displacement")
		# Solve variational problem
		solve(Form_s == 0, self.d, self.bcs, J = Gain_s, \
			solver_parameters  = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
		end()

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"
