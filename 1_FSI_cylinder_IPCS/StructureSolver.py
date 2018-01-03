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

		self.V_element = VectorElement(self.ElementType, self.mesh.ufl_cell() , self.ElementDegree)

		# Mixed space for computing both displacement and rate of change of displacement
		self.V_element = VectorElement(self.ElementType, self.mesh.ufl_cell() , self.ElementDegree)
		#self.V2_temp = self.V_element*self.V_element

		#self.V2_space = FunctionSpace(self.mesh, self.V2_temp)
		self.V2_space = FunctionSpace(self.mesh, MixedElement([self.V_element, self.V_element]))
		self.V2_space = FunctionSpace(self.mesh, self.V_element*self.V_element)

		self.V = TestFunction(self.V2_space)
		self.dU = TrialFunction(self.V2_space)
		self.U = Function(self.V2_space)
		self.U0 = Function(self.V2_space)

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

		## Functions for solver
		#self.d_ = Function(self.V_space) 		# displacement from current iteration
		#self.d_n = Function(self.V_space)		# displacement from prevous iteration
		#self.v_ = Function(self.V_space)		# velocity from the current iteration
		#self.v_n = Function(self.V_space)		# velocity from the previous iteration

		## do I need these for velocity too?
		#self.du = TrialFunction(self.V_space)	# incemental displacement
		#self.v = TestFunction(self.V_space)		# Test function



###########################
	#	self.d00_s = Function(self.V_space)
###########################
	def Structure_Problem_Solver(self, p_s, F):

		if self.solver == "Linear":
			self.Linear_Elastic_Solver(p_s, F)
		elif self.solver == "NeoHookean":
			#self.Incompressible_NeoHookean_Solver(p_s, F)
			self.Compressible_NeoHookean_Solver(p_s, F)
		elif self.solver =="St_venant":
			self.Compressible_St_venant(p_s, F)

		else:
			print "Error. The only solvers available for the structure are Linear or NeoHookean"

		#self.d_dot = project( (self.d_ - self.d_n)/p_s.dt, self.V_space, solver_type = "mumps", \
			#form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

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
		#self.T_hat = Constant((0.0,0.0))

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

	def Compressible_St_venant(self, p_s, F):

		print ""
		print ""
		print "ENTERING STRUCTURE St Venant SOLVER"
		print ''
		print ''

		# displacements and velocities at mid points
		self.u_mid = 0.5*(self.u0 + self.u)
		self.v_mid = 0.5*(self.v0 + self.v)

		I = Identity(p_s.Dim)

		# Project stress from fluid tensor space onto structure tensor space
		self.sigma_FSI = project(F.sigma_FSI, self.T_space, solver_type = "mumps",\
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		#project seems wrong when values on boundary are printed...
		#I = Identity(self.d_.cell().d)

		# Kinematics
		# Only use 2nd pk stres evaluated at u_mid it seems.
		self.F =  variable(I + grad(self.u_mid))		# Deformation gradient tensor
		self.C =  variable(self.F.T*self.F)			# Right Cauchy-Green tensor

		# Green lagrange stress tensor (green-st venant strain tensor)
		self.E = variable(0.5*(self.C-I))

		# stored strain energy
		self.psi = p_s.lambda_s/2*(tr(self.E)**2)+p_s.mu_s*tr(self.E*self.E)

		# 2nd piola kirchoff stress tensor
		self.S2 = diff(self.psi, self.E)

		# 1st piola kirchoff stress tensor
		#evaluated at u_mid
		self.S1 = self.F*self.S2

		# Invariants of deformation tensor
		self.J = det(self.F)

		# Evaluate

		self.k = Constant(p_s.dt)
		# Variational form for hyperelasticity

		# Traction on boundary

		# solves and gives 0 displacment
		#self.T_hat = Constant((0.0, 0.0))

		# Solves and gives too great a displacement.
		n_f = FacetNormal(F.mesh)
		n_s = FacetNormal(self.mesh)

		# To be consistent, use u_mid for T_hat too...
		self.T_hat = (self.J*inv(self.F)*self.sigma_FSI).T * n_s
		# gives same results as index method
		#self.T_hat = Constant((0.0, 0.0))
		#self.T_hat = as_tensor( self.J*inv(self.F)[k, j]*self.sigma_FSI[j, i]*n_s[k] , (i, ) )

		#self.Dim = self.mesh.topology().dim()
		#i, j, k, l, m = indices(5)
		#self.delta = Identity(self.Dim)


		# Piola map
		#if p_s.t <= 10*p_s.dt:
		#	self.T_hat = Constant((0.0, 0.0))

		# solves and gives 0 displacment
		#self.T_hat = Constant((0.0, 0.0))

		# Total potential energy
		#self.Pi = self.psi*self.dx - dot(self.T_hat, self.d_)*self.ds(2) - dot(self.B, self.d_)*self.dx

		# The variational form corresponding to hyperelasticity

		self.L = p_s.rho_s*inner(self.v - self.v0, self.u_t)*self.dx \
		+ self.k*inner(self.S1, grad(self.u_t))*self.dx \
		- self.k*inner(self.B, self.u_t)*self.dx \
		+ inner(self.u - self.u0, self.v_t)*self.dx \
		- self.k*inner(self.v_mid, self.v_t)*self.dx

		self.L = self.L-self.k*inner(self.T_hat, self.u_t)*self.ds(2)
		#begin("Computing structure displacement")
		# get Neumann BCs on the stress. This is the hard part.
		#Neumann stuff
		# in problem_definitions:
		#neumann_conditions - return neumann boundary conditions for the stress field
		#neumann_boundaries - return boundaries over which Neumann conditions act
		# T_hat expression... try with constant for now.
		#self.T_hat = Constant((0.0, 0.0))
		#self.T_hat = Expression(('0.0', '0.0'))

		# In Cylinder the boundary S.fsi is defined and marked as 2.
		# the definition here of ds differes a little from their one.

		#compiled_boundary = compile_subdomains(S.fsi)
		#compiled_boundary.mark(boundary,7)


		# introduced in line 499 of solution_algorithms
		self.a = derivative(self.L, self.U, self.dU)
		# Setup problem
		# may have to fix bcu to align with Twist.
		# a bit iffy on how solve, step, update interact in CBC twist... write out, comb over.

		problem = NonlinearVariationalProblem(self.L, self.U, self.bcs, self.a)
		solver = NonlinearVariationalSolver(problem)
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["maximum_iterations"] = 100

		solver.solve()


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

		#self.T_hat = Constant((0.0,0.0))

		self.n = FacetNormal(self.mesh)
		#self.T_hat = dot(self.sigma_FSI, self.n)

		# Total potential energy
		#self.Pi = self.psi*self.dx - dot(self.T_hat, self.d_)*self.ds(2) - dot(self.B, self.d_)*self.dx
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
