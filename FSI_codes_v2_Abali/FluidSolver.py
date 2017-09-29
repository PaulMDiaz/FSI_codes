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
		self.T_space = TensorFunctionSpace(self.mesh, self.ElementType, self.VelocityElementDegree)

		## Functions for results
		self.v_ = Function(self.V_space, name = 'u')
		self.p_ = Function(self.S_space, name = 'p')
                # called only for saving results... weird...

		# Functions for solver at previos and current time steps... could also be results...
		self.u1 = Function(self.V_space)
		self.p1 = Function(self.S_space)
		self.u0 = Function(self.V_space)
		self.u_mesh = Function(self.V_space)
		#if self.solver == "Chorin":
			#self.us = Function(self.V_space)	#not sure what this is for?
		#else:
		#	print "Error. The solver for the fluid problem should be set to Chorin"

		self.sigma_FSI = Function(self.T_space) # pressure on interface. a function..

		# Trial and test functions
		self.u = TrialFunction(self.V_space)
		self.p = TrialFunction(self.S_space)
		self.du = TestFunction(self.V_space)
		self.dp = TestFunction(self.S_space)

		# Trial and test functions for mesh motion
		self.del_u_mesh = TestFunction(self.V_space)
		self.du_mesh = TrialFunction(self.V_space)


	def Fluid_Problem_Solver(self, DC, S):

		if self.solver == "Chorin":
			self.Chorin_Fluid_Solver(DC, S)

		else:
			print "Error. The solver for the fluid problem should be set to Chorin"

	def Chorin_Fluid_Solver(self, DC, S):
		print ""
		print ""
		print "ENTERING FLUID SOLVER"
		print ''
		print ''

		I = Identity(DC.Dim)

                # Verify variational form

		# Tentative velocity step
                # should u mesh appear anywhere else in this??
		# In this case u represents the tentative velocity

		 # temporal term. rate of change of velocity
		 # convective term (taken out to match paper?)
		 # diffusive term.

		 # Variational form for step 1: Tentative velocity
		#F1 = (1/DC.dt)*inner(self.u - self.u0, self.du)*self.dx\
                    # + inner(grad(self.u0)*(self.u0 - self.u_mesh), self.du)*self.dx\
                     #+ DC.nu_f*inner(grad(self.u), grad(self.du))*self.dx\
                     #- inner(self.f, self.du)*self.dx
		# F1 variational form matches demo except for u_mesh term.
		# a2 and a3 match perfectly
		F1 = (1/DC.dt)*inner(self.u - self.u0, self.du)*self.dx\
                      + inner(grad(self.u0)*(self.u0-self.u_mesh), self.du)*self.dx\
                      + DC.nu_f*inner(grad(self.u), grad(self.du))*self.dx\
                      - inner(self.f, self.du)*self.dx
		self.a1 = lhs(F1)
		self.L1 = rhs(F1)

		# Step 2: Pressure update
		self.a2 = inner(grad(self.p), grad(self.dp))*self.dx
		self.L2 = -(1/DC.dt)*div(self.u1)*self.dp*self.dx #u1 was us... defined with chorin

		# Step 3: Velocity update
		self.a3 = inner(self.u, self.du)*self.dx
		self.L3 = inner(self.u1, self.du)*self.dx - DC.dt*inner(grad(self.p1), self.du)*self.dx #u1 was us... defined with chorin

		# Assemble matrices
		self.A1 = assemble(self.a1)
		self.A2 = assemble(self.a2)
		self.A3 = assemble(self.a3)

		# Define the matrices to solve the system

		# Apply boundary conditions to matrices
		##[bc.apply(self.A1) for bc in self.bcu]
		##[bc.apply(self.A2) for bc in self.bcp]

		# Compute tentative velocity step
		begin("Computing tentative velocity")
		self.b1 = assemble(self.L1)
		[bc.apply(self.A1, self.b1) for bc in self.bcu]
		##self.b1 = assemble(self.L1)
		##[bc.apply(self.b1) for bc in self.bcu]
		##solve(self.A1, self.u1.vector(), self.b1, "gmres", "default")
		solve(self.A1, self.u1.vector(), self.b1, "gmres", "ilu")
		end()

		prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
		# parameters['krylov_solver']['nonzero_initial_guess'] = True

		# Pressure correction
		begin("Computing pressure correction")
		self.b2 = assemble(self.L2)
		##self.b2 = assemble(self.L2)
		##[bc.apply(self.b2) for bc in self.bcp]
		#solve(self.A2, self.p1.vector(), self.b2, "gmres", "default")
		solve(self.A2, self.p1.vector(), self.b2, "gmres", "ilu")
		end()

		# Velocity correction
		begin("Computing velocity correction")
		self.b3 = assemble(self.L3)
		[bc.apply(self.A3, self.b3) for bc in self.bcu]
		##self.b3 = assemble(self.L3)
		##solve(self.A3, self.u1.vector(), self.b3, "gmres", "default")
		solve(self.A3, self.u1.vector(), self.b3, "gmres", "ilu")
		end()

		self.tau = DC.mu_f*(grad(self.u1) + grad(self.u1).T)

        # sigma_FSI stress at boundary
		self.sigma_FSI = project(-self.p1*I + self.tau, self.T_space, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		print ""
		print ""
		print "EXITING FLUID SOLVER"
