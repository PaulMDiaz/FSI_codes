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
###########################
		self.d00_s = Function(self.V_space)
###########################
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

		# project fluid stress onto structure tensor space.
		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		I = Identity(DC.Dim)

		self.sigma = 2.0*DC.mu_s*sym(grad(self.u)) + DC.lambda_s*tr(sym(grad(self.u)))*I

		self.F = I + grad(self.u)
		self.J = det(self.F)

		#n = FacetNormal(self.mesh)

		self.T_fat = self.J*inv(self.F)*self.sigma_FSI*self.N
		#self.T_fat = dot(self.J*inv(self.F)*self.sigma_FSI,self.N)
		#T_fat = dot(S.J*inv(S.F)*S.sigma_FSI,n)

		#self.T_hat = Constant((0.0,0.0))
		n = FacetNormal(self.mesh)
		self.T_hat = dot(self.sigma_FSI, n)

		print ""
		print ""
		print "ENTERING STRUCTURE LINEAR SOLVER"
		print ''
		print ''

		# Setup variational problem

		# dA is ds, dV is dx
		# not certain about this a line differs from demo.
		self.a = inner(self.sigma, grad(self.du))*self.dV
		self.L = inner(self.T_hat, self.du)*self.dA(3) # dA(3) = ds(3) = FSI boundary

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

		self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N

		#self.n = FacetNormal(self.mesh)
		#self.T_hat = dot(self.sigma_FSI, self.n)
		#self.T_fat = dot(dot(self.J*inv(self.F),self.sigma_FSI), self.n)

		#n = FacetNormal(self.mesh)
		#self.T_hat = dot(self.sigma_FSI, n)

		#self.T_hat = Constant((0.0,0.0))

		# Total potential energy
		# dA(3) integrate on interface only or perhaps wrong form...
		#self.Pi = self.psi*self.dV - dot(self.T_hat, self.d)*self.dA(3) - dot(self.B, self.d)*self.dV
		self.Pi = self.psi*self.dV - dot(self.T_hat, self.d)*self.dA(3) - dot(self.B, self.d)*self.dV

		# First directional derivative of Pi about d in the direction of v
		Form_s = derivative(self.Pi, self.d, self.du)




		#######################################################################
		# Alternative using code from Abali. Don't entirely understand. indices were not used in fluid solver or other parts of this code...
		# requires edits to other parts too in order to update with time.
		# using indices like Abali's code:
		self.n = FacetNormal(self.mesh)

		self.Dim = self.mesh.topology().dim()
		i,j, k, l, m = indices(5)
		self.delta = Identity(self.Dim)

		self.F_s = as_tensor( self.d[k].dx(i) + self.delta[k, i], (k, i) )
		self.J_s = det(self.F_s)

		self.C_s = as_tensor( self.F_s[k, i]*self.F_s[k, j], (i, j) )
		self.E_s = as_tensor(1./2.*(self.C_s[i, j] - self.delta[i, j]), (i, j) )
		self.S_s = as_tensor( DC.lambda_s*self.E_s[k, k]*self.delta[i, j] + 2.*DC.mu_s*self.E_s[i, j], (i, j))

		# S_s = as_tensor( lambda_s*E[k, k]*delta[i, j] + 2.*mu_s*E[i, j], (i, j))
		self.P_s = as_tensor( self.F_s[i, j]*self.S_s[k, j], (k, i) )

		self.t_hat = as_tensor( self.J_s*inv(self.F_s)[k, j]*self.sigma_FSI[j, i]*self.n[k] , (i, ) )


		Form_s = ( DC.rho_s*(self.d-2.*self.d0+self.d00_s)[i]/(DC.dt*DC.dt)*self.du[i] + self.P_s[k, i]*self.du[i].dx(k) - DC.rho_s*self.B[i]*self.du[i] )*self.dV - \
		     self.t_hat[i]*self.du[i]*self.dA(3)

		########################################################################


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
		#self.psi = (DC.mu_s/2)*(self.Ic - 3) - DC.mu_s*ln(self.J) + (DC.lambda_s/2)*(ln(self.J))**2
		self.psi = (DC.mu_s/2)*(self.Ic-3)

		self.T_hat = self.J*inv(self.F)*self.sigma_FSI*self.N
		#self.T_hat = Constant((0.0,0.0))

		#self.n = FacetNormal(self.mesh)
		#self.T_hat = dot(self.sigma_FSI, self.n)

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
