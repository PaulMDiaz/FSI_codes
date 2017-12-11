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
	def __init__(self, Mesh, ElementType, VelocityElementDegree, PressureElementDegree, FluidSolverMethod, FluidBodyForce):

		# Initialize mesh and solver parameters
		self.mesh = Mesh
		self.ElementType = ElementType
		self.VelocityElementDegree = VelocityElementDegree
		self.PressureElementDegree = PressureElementDegree
		self.solver = FluidSolverMethod
		self.f = FluidBodyForce

		# Variational spaces
		self.S_space = FunctionSpace(self.mesh, self.ElementType, self.PressureElementDegree)
		self.V_space = VectorFunctionSpace(self.mesh, self.ElementType, self.VelocityElementDegree)
		self.T_space = TensorFunctionSpace(self.mesh, self.ElementType, self.PressureElementDegree)

		## Functions for results
		self.u_res = Function(self.V_space, name = 'u')
		self.p_res = Function(self.S_space, name = 'p')
                # called only for saving results... weird...

		# Functions for solver at previos and current time steps... could also be results...
		self.u_n = Function(self.V_space) 	# Previous time step
		self.u_  = Function(self.V_space)	# Current time step
		self.p_n = Function(self.S_space)	# Previous time step
		self.p_  = Function(self.S_space)	# Current time step

		# load data
		#u_temp = np.loadtxt('nodal_u_64')
		#self.u_n.vector()[:] = u_temp
		#p_temp = np.loadtxt('nodal_p_64')
		#self.p_n.vector()[:] = p_temp


		# check these
		self.u0 = Function(self.V_space)	# doesn't appear to be used anywhere... check
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


	def Fluid_Problem_Solver(self, p_s, S):

		if self.solver == "Chorin":
			self.Chorin_Fluid_Solver(p_s, S)

		else:
			print "Error. The solver for the fluid problem should be set to Chorin"

	def Chorin_Fluid_Solver(self, p_s, S):
		print ""
		print ""
		print "ENTERING FLUID SOLVER"
		print ''
		print ''

		I = Identity(p_s.Dim)

		print(p_s.rho_f)
		print(p_s.mu_f)
		print(p_s.nu_f)

		# Define expressions used in variational forms
		self.U  = 0.5*(self.u_n + self.u)
		self.n  = FacetNormal(self.mesh)
		self.f  = Constant((0, 0))
		self.k  = Constant(p_s.dt)
		self.mu_f = Constant(p_s.mu_f)
		self.rho_f = Constant(p_s.rho_f)


		# Define symmetric gradient
		def epsilon(u):
		    return sym(nabla_grad(u))

		# Define stress tensor
		def sigma(u, p):
		    return 2*self.mu_f*epsilon(u) - p*Identity(len(u))

		self.F1 = self.rho_f*dot((self.u - self.u_n) / self.k, self.v)*self.dx \
		   + self.rho_f*dot(dot((self.u_n-self.u_mesh), nabla_grad(self.u_n)), self.v)*self.dx \
		   + inner(sigma(self.U, self.p_n), epsilon(self.v))*self.dx \
		   + dot(self.p_n*self.n, self.v)*self.ds \
		   - dot(self.mu_f*nabla_grad(self.U)*self.n, self.v)*self.ds \
		   - dot(self.f, self.v)*self.dx
		#self.F1 = (1/self.k)*inner(self.u - self.u_n, self.v)*self.dx\
        #             + inner(grad(self.u_n)*(self.u_n - self.u_mesh), self.v)*self.dx\
        #             + p_s.nu_f*inner(grad(self.u), grad(self.v))*self.dx\
        #             - inner(self.f, self.v)*self.dx
		self.a1 = lhs(self.F1)
		self.L1 = rhs(self.F1)

		# Pressure update
		self.a2 = dot(nabla_grad(self.p), nabla_grad(self.q))*self.dx
		self.L2 = dot(nabla_grad(self.p_n), nabla_grad(self.q))*self.dx - (1/self.k)*div(self.u_)*self.q*self.dx

		# Velocity update
		self.a3 = dot(self.u, self.v)*self.dx
		self.L3 = dot(self.u_, self.v)*self.dx - self.k*dot(nabla_grad(self.p_ - self.p_n), self.v)*self.dx

		# Assemble matrices
		self.A1 = assemble(self.a1)
		self.A2 = assemble(self.a2)
		self.A3 = assemble(self.a3)

		# Define the matrices to solve the system

		[bc.apply(self.A1) for bc in self.bcu]
		[bc.apply(self.A2) for bc in self.bcp]
		# maybe this line?
		#[bc.apply(self.A3) for bc in self.bcu]

		# Compute tentative velocity step
		#'hypre_amg'
		begin("Computing tentative velocity")
		self.b1 = assemble(self.L1)
		#[bc.apply(self.b1) for bc in self.bcu]
		[bc.apply(self.A1,self.b1) for bc in self.bcu]
		solve(self.A1, self.u_.vector(), self.b1, 'bicgstab', 'ilu')
		end()
		#print('u test 1:',self.u_(0.1,0.2))
		#prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
		# parameters['krylov_solver']['nonzero_initial_guess'] = True

		# Pressure correction
		#'hypre_amg'
		begin("Computing pressure correction")
		self.b2 = assemble(self.L2)
		[bc.apply(self.b2) for bc in self.bcp]
		solve(self.A2, self.p_.vector(), self.b2, 'bicgstab', 'ilu')
		end()
		#print('p test 1:',self.p_(0.15,0.2))
		#'hypre_amg' or 'imu'
		# Velocity correction
		begin("Computing velocity correction")
		#, 'sor'
		self.b3 = assemble(self.L3)
		# this line is not in earlier code... seems to assert BC.
		#[bc.apply(self.A3,self.b3) for bc in self.bcu]
		solve(self.A3, self.u_.vector(), self.b3, 'cg')
		end()
		#print('u test 2:',self.u_(0.1,0.2))
		self.tau = p_s.mu_f*(grad(self.u_) + grad(self.u_).T)

        # sigma... seems unnecessaryily cumbersome to compute stress over entire domain.
		self.sigma_FSI = project(-self.p_*I + self.tau, self.T_space, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		print ""
		print ""
		print "EXITING FLUID SOLVER"
