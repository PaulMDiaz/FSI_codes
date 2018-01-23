#  Configures variational problem for incompressible Navier
#  Stokes equations to solve for pressure and velocity distributions
#  in a Driven Cavity Flow.
#
#  Created: 12 May 2015
#  Modified: 5 June 2015
#
#  Inputs:
#    IOInfo - An object which contains all inputs and outputs for
#      the Navier Stokes solver implemented herein

#from dolfin import *
from fenics import *
import numpy as np
import pylab as plt

class Fluid_Solver(object):
	"""docstring for ClassName"""
	def __init__(self, Mesh, ElementType, VelocityElementDegree, PressureElementDegree, FluidSolverMethod, FluidBodyForce, p_s):

		# Initialize mesh and solver parameters
		self.mesh = Mesh
		self.ElementType = ElementType
		self.VelocityElementDegree = VelocityElementDegree
		self.PressureElementDegree = PressureElementDegree
		self.solverMethod = FluidSolverMethod
		self.f = FluidBodyForce

		# Variational spaces
		self.S_space = FunctionSpace(self.mesh, self.ElementType, self.PressureElementDegree)
		self.V_space = VectorFunctionSpace(self.mesh, self.ElementType, self.VelocityElementDegree)
		# Abali has all degrees as 1. Doesn't seem quite right.
		#self.V_space = VectorFunctionSpace(self.mesh, self.ElementType, self.PressureElementDegree)
		self.T_space = TensorFunctionSpace(self.mesh, self.ElementType, self.PressureElementDegree)

		## Functions for results
		self.u_res = Function(self.V_space, name = 'u')
		self.p_res = Function(self.S_space, name = 'p')
        # called only for saving results... weird...

		# Functions for solver at previos and current time steps... could also be results...
		self.u0 = Function(self.V_space) 	# Previous time step
		self.u1  = Function(self.V_space)	# Current time step
		self.p1 = Function(self.S_space)	# Current time step
		self.p0  = Function(self.S_space)	# Previous time step

		# load data
		#u_temp = np.loadtxt('nodal_u_64')
		#self.u0.vector()[:] = u_temp
		#p_temp = np.loadtxt('nodal_p_64')
		#self.p0.vector()[:] = p_temp


		# check these
	#	self.u0 = Function(self.V_space)	# doesn't appear to be used anywhere... check
		self.u_mesh = Function(self.V_space)# Current time step

		# self.sigma_FSI = Function(self.T_space) # pressure on interface. a function..

		# Trial and test functions
		# Define trial and test functions
		self.u = TrialFunction(self.V_space)
		self.p = TrialFunction(self.S_space)
		self.v = TestFunction(self.V_space)
		self.q = TestFunction(self.S_space)

		# Trial and test functions for mesh motion
		self.del_u_mesh = TestFunction(self.V_space)
		self.du_mesh = TrialFunction(self.V_space)

		# Define coefficients
		self.k = Constant(p_s.dt)
		self.U  = 0.5*(self.u0 + self.u)
		self.n  = FacetNormal(self.mesh)

		self.cells = CellFunction("size_t", self.mesh)
		self.facets = FacetFunction("size_t", self.mesh)
		self.dx = Measure('dx', domain = self.mesh, subdomain_data = self.cells)
		self.ds = Measure('ds', domain = self.mesh, subdomain_data = self.facets)

		# beta is a flag. set to zero when periodic boundary conditions are used.
		self.beta = 1

		# Define symmetric gradient
		def epsilon(u):
		    return sym(nabla_grad(u))

		# Define stress tensor
		def sigma(u, p):
		    return 2*p_s.mu_f*epsilon(u) - p*Identity(len(u))

		# Think that u_mesh subtraction is in correct place. Not inside grad.
		# Tentative velocity step
		self.F1 = p_s.rho_f*(1/self.k)*inner(self.u - self.u0, self.v)*self.dx\
			+ p_s.rho_f*inner(grad(self.u0)*(self.u0-self.u_mesh), self.v)*self.dx\
			+ inner(sigma(self.U, self.p0), epsilon(self.v))*self.dx \
			+ inner(self.p0*self.n, self.v)*self.ds \
			- self.beta*p_s.nu_f*inner(grad(self.U).T*self.n, self.v)*self.ds \
			- inner(self.f, self.v)*self.dx
		self.a1 = lhs(self.F1)
		self.L1 = rhs(self.F1)

		# Pressure correction
		self.a2 = inner(grad(self.p), grad(self.q))*self.dx
		self.L2 = inner(grad(self.p0), grad(self.q))*self.dx \
			-(1.0/self.k)*div(self.u1)*self.q*self.dx
			#-div(self.u1)*self.q*self.dx
			#-(1.0/self.k)*div(self.u1)*self.q*self.dx

		# Velocity correction
		self.a3 = inner(self.u, self.v)*self.dx
		self.L3 = inner(self.u1, self.v)*self.dx\
		     - self.k*inner(grad(self.p1-self.p0), self.v)*self.dx

		# Assemble matrices
		self.A1 = assemble(self.a1)
		self.A2 = assemble(self.a2)
		self.A3 = assemble(self.a3)

		# Use amg preconditioner if available
		self.prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

		# Use nonzero guesses - essential for CG with non-symmetric BC
		parameters["std_out_all_processes"] = False;
		parameters['krylov_solver']['nonzero_initial_guess'] = True

		# Create files for storing solution
		#ufile = File("results_cylinder/velocity.pvd")
		#pfile = File("results_cylinder/pressure.pvd")

		[bc.apply(self.A1) for bc in self.bcu]
		[bc.apply(self.A2) for bc in self.bcp]

	def solve_fluid_problem(self, p_s):
		if self.solverMethod == "IPCS":
			self.IPCS_Fluid_Solver(p_s)

		else:
			print "Error. The solver for the fluid problem should be set to IPCS"

	def IPCS_Fluid_Solver(self, p_s):
		print ""
		print ""
		print "ENTERING FLUID SOLVER"
		print ''
		print ''

		I = Identity(p_s.Dim)

		# Compute tentative velocity step
		begin("Computing tentative velocity")
		self.b1 = assemble(self.L1)
		#[bc.apply(self.A1, self.b1) for bc in self.bcu]
		[bc.apply(self.b1) for bc in self.bcu]
		solve(self.A1, self.u1.vector(), self.b1, "bicgstab", "default")
		end()

		# Pressure correction
		begin("Computing pressure correction")
		self.b2 = assemble(self.L2)
		[bc.apply(self.b2) for bc in self.bcp]
		solve(self.A2, self.p1.vector(), self.b2, "bicgstab", self.prec)
		end()

		# Velocity correction
		begin("Computing velocity correction")
		self.b3 = assemble(self.L3)
		# Navier stokes cylinder demo does not enforce BC here. Not condident on conrrectedness of that approach. 
		#[bc.apply(self.A3, self.b3) for bc in self.bcu]
		solve(self.A3, self.u1.vector(), self.b3, "bicgstab", "default")

		#print('u test 2:',self.u1(0.1,0.2))
		self.tau = p_s.mu_f*(grad(self.u1) + grad(self.u1).T)

        # Sigma
		self.sigma_FSI = project(-self.p1*I + self.tau, self.T_space, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )
		#self.sigma_FSI = project(-self.p1*I + self.tau, self.T_space, solver_type = "mumps", \
			#form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
		print ""
		print ""
		print "EXITING FLUID SOLVER"
