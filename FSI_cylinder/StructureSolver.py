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


		## Function for results
		self.d_res = Function(self.V_space, name = 'd')

		# Functions for solver
		self.d_ = Function(self.V_space) 		# displacement from current iteration
		self.d_n = Function(self.V_space)		# displacement from prevous iteration
		self.du = TrialFunction(self.V_space)	# incemental displacement
		self.v = TestFunction(self.V_space)		# Test function

###########################
	#	self.d00_s = Function(self.V_space)
###########################
	def Structure_Problem_Solver(self, p_s, F):

		if self.solver == "Linear":
			self.Linear_Elastic_Solver(p_s, F)
		elif self.solver == "NeoHookean":
			#self.Incompressible_NeoHookean_Solver(p_s, F)
			self.Compressible_NeoHookean_Solver(p_s, F)
		else:
			print "Error. The only solvers available for the structure are Linear or NeoHookean"

		self.d_dot = project( (self.d_ - self.d_n)/p_s.dt, self.V_space, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

	def Linear_Elastic_Solver(self, p_s, F):

		# project fluid stress onto structure tensor space.
		# don't know that this succeeds.

		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		I = Identity(p_s.Dim)

		self.sigma = 2.0*p_s.mu_s*sym(grad(self.du)) + p_s.lambda_s*tr(sym(grad(self.du)))*I

		self.F = I + grad(self.du)
		self.J = det(self.F)

		#n = FacetNormal(self.mesh)

		#self.T_fat = self.J*inv(self.F)*self.sigma_FSI*self.N
		#self.T_fat = dot(self.J*inv(self.F)*self.sigma_FSI,self.N)
		#T_fat = dot(S.J*inv(S.F)*S.sigma_FSI,n)

		#self.T_hat = Constant((0.0,0.0))
		n = FacetNormal(self.mesh)
		self.T_hat = dot(self.sigma_FSI, n)
		self.T_hat = Constant((0.0,0.0))

		print ""
		print ""
		print "ENTERING STRUCTURE LINEAR SOLVER"
		print ''
		print ''

		# Setup variational problem

		# ds is ds, dx is dx
		# not certain about this a line differs from demo.
		self.a = inner(self.sigma, grad(self.v))*self.dx
		self.L = inner(self.T_hat, self.v)*self.ds(2) # ds(3) = ds(3) = FSI boundary

		# Solve variational problem
		begin("Computing structure displacement")
		solve(self.a == self.L, self.d_, self.bcs, solver_parameters={"symmetric":True})
		end()
			#self.Incompressible_NeoHookean_Solver(p_s, F)

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"

	def Compressible_NeoHookean_Solver(self, p_s, F):

		print ""
		print ""
		print "ENTERING STRUCTURE NEOHOOKEAN SOLVER"
		print ''
		print ''

		I = Identity(p_s.Dim)

		# Project stress from fluid tensor space onto structure tensor space
		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		#project seems wrong when values on boundary are printed...

		# Kinematics
		self.F = I + grad(self.d_)			# Deformation gradient tensor
		self.C = self.F.T*self.F			# Right Cauchy-Green tensor

		# Invariants of deformation tensor
		self.Ic = tr(self.C)
		self.J = det(self.F)

		# Elasticity parameters E and nu should be defined

		# Stored strain energy density (compressible neo-Hookean model)
		self.psi = (p_s.mu_s/2)*(self.Ic - 3) - p_s.mu_s*ln(self.J) + (p_s.lambda_s/2)*(ln(self.J))**2

		#original line
		#self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N

		# Traction force (per unit reference area)
		# (first piola kirchoff stress tensor, also called called lagrangian stress tensor)
		self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.n

		#self.T_fat = self.J*inv(self.F)*self.sigma_FSI*self.N
		#self.n_f = FacetNormal(F.mesh)

		#self.T_hat = -self.J*inv(self.F)*F.sigma_FSI*self.N
		#self.T_fat = -dot(self.J*F.sigma_FSI*inv(self.F).T,self.n_f)
		#self.T_fat = -self.J*F.sigma_FSI*inv(self.F).T*self.n_f


		#self.T_fat_2 = project(self.T_fat, self.V_space, solver_type = "mumps",\
		form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2}

		#self.n = FaceTtNormal(self.mesh)
		#self.T_hat = dot(self.sigma_FSI, self.n)
		#self.T_fat = dot(dot(self.J*inv(self.F),self.sigma_FSI), self.n)

		#n = FacetNormal(self.mesh)
		#self.T_hat = dot(self.sigma_FSI, n)

		#self.T_hat = Constant((0.0,1.0))

		# Total potential energy
		# ds(3) integrate on interface only or perhaps wrong form...
		#self.Pi = self.psi*self.dx - dot(self.T_hat, self.d)*self.ds(3) - dot(self.B, self.d)*self.dx
		self.Pi = self.psi*self.dx - dot(self.T_hat, self.d_)*self.ds(2) - dot(self.B, self.d_)*self.dx

		# First directional derivative of Pi about d in the direction of v
		Form_s = derivative(self.Pi, self.d_, self.v)

		# Jacobian of the directional derivative Fd
		Gain_s = derivative(Form_s, self.d_, self.du)

		begin("Computing structure displacement")
		# Solve variational problem. Why these solver parameters?

		#Could use:
		parameters["form_compiler"]["cpp_optimize"] = True
		ffc_options = {"optimize": True}

		#solve(Form_s == 0, self.d, self.bcs, J=Gain_s,
		#      form_compiler_parameters=ffc_options)

		solve(Form_s == 0, self.d_, self.bcs, J = Gain_s, \
			solver_parameters  = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
		end()

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"

	def Incompressible_NeoHookean_Solver(self, p_s, F):

		print ""
		print ""
		print "ENTERING STRUCTURE NEOHOOKEAN SOLVER"
		print ''
		print ''

		I = Identity(p_s.Dim)

		# Project stress from fluid tensor space onto structure tensor space.
		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		# Kinematics
		self.F = I + grad(self.d_)			# Deformation gradient
		self.C = self.F.T*self.F			# Right Cauchy-Green tensor

		# Invariants of deformation tensor
		self.Ic = tr(self.C)
		self.J = det(self.F)

		# Elasticity parameters E and nu should be defined.

		# Stored strain energy density
		#self.psi = (p_s.mu_s/2)*(self.Ic - 3) - p_s.mu_s*ln(self.J) + (p_s.lambda_s/2)*(ln(self.J))**2
		self.psi = (p_s.mu_s/2)*(self.Ic-3)

		#self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N
		self.T_hat = Constant((0.0,0.0))

		#self.n = FacetNormal(self.mesh)
		#self.T_hat = dot(self.sigma_FSI, self.n)

		# Total potential energy
		self.Pi = self.psi*self.dx - dot(self.T_hat, self.d_)*self.ds(2) - dot(self.B, self.d_)*self.dx

		# First directional derivative of Pi about d in the direction of v
		Form_s = derivative(self.Pi, self.d_, self.v)

		# Jacobian of the directional derivative Fd
		Gain_s = derivative(Form_s, self.d_, self.du)

		begin("Computing structure displacement")
		# Solve variational problem
		solve(Form_s == 0, self.d_, self.bcs, J = Gain_s, \
			solver_parameters  = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
		end()

		print ""
		print ""
		print "EXITING STRUCTURE SOLVER"
