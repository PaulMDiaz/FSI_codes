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

class FluidSolver(object):
	"""docstring for ClassName"""
	def __init__(self, mesh, elementType, velocityElementDegree, pressureElementDegree, fluidSolverMethod, fluidBodyForce, p_s, meshSolver):

		# Initialize mesh and solver parameters
		self.mesh = mesh
		self.elementType = elementType
		self.velocityElementDegree = velocityElementDegree
		self.pressureElementDegree = pressureElementDegree
		self.solverMethod = fluidSolverMethod
		self.f = fluidBodyForce

		# Variational spaces
		self.scalarSpace = FunctionSpace(self.mesh, self.elementType, self.pressureElementDegree)
		self.vectorSpace = VectorFunctionSpace(self.mesh, self.elementType, self.velocityElementDegree)
		self.tensorSpace = TensorFunctionSpace(self.mesh, self.elementType, self.pressureElementDegree)
		self.tensorSpaceRef = TensorFunctionSpace(p_s.meshRef, self.elementType, self.pressureElementDegree)

		## Functions for results
		self.u_res = Function(self.vectorSpace, name = 'u')
		self.p_res = Function(self.scalarSpace, name = 'p')
        # called only for saving results... weird...

		# Functions for solver at previos and current time steps... could also be results...
		self.u0 = Function(self.vectorSpace) 	# Previous time step
		self.u1  = Function(self.vectorSpace)	# Current time step
		self.p1 = Function(self.scalarSpace)	# Current time step
		self.p0  = Function(self.scalarSpace)	# Previous time step

		# Load initial conditions to u0 and v0. Otherwise set to 0.
		#self.u0 = Constant((0,)*self.vectorSpace.mesh().geometry().dim())
		#self.v0 = Constant((0,)*self.vectorSpace.mesh().geometry().dim())
		#self.u0 = Constant((0,)*2)
		#self.v0 = Constant((0,)*2)

		# load data
		#u_temp = np.loadtxt('nodal_u_64')
		#self.u0.vector()[:] = u_temp
		#p_temp = np.loadtxt('nodal_p_64')
		#self.p0.vector()[:] = p_temp

		# fluid velocity on interface, functionspace and function
		self.interfaceFluidVectorSpace = VectorFunctionSpace(p_s.interfaceFluid, "CG", 1)
		interfaceFluidCoords = self.interfaceFluidVectorSpace.tabulate_dof_coordinates().reshape((self.interfaceFluidVectorSpace.dim(),-1))

		# Function for fluid interface velocity
		self.fluidInterfaceVelocity = Function(self.interfaceFluidVectorSpace)
		self.fluidInterfaceVelocity.vector().zero()

		# Vector function space degree 1 for copying mesh velocity to:
		self.vectorSpaceMesh = VectorFunctionSpace(self.mesh, self.elementType, 1)
		self.meshVelocityCurrent1 = Function(self.vectorSpaceMesh) 	# Previous time step

		# self.interfaceFluidTensorSpace = TensorFunctionSpace(p_s.interfaceFluid, "CG", 1)
		# interfaceFluidTensorCoords = self.interfaceFluidTensorSpace.tabulate_dof_coordinates().reshape((self.interfaceFluidTensorSpace.dim(),-1))

		# self.interfaceRefTensorSpace = TensorFunctionSpace(p_s.interfaceMesh, "CG", 1)
		# interfaceRefTensorCoords = self.interfaceRefTensorSpace.tabulate_dof_coordinates().reshape((self.interfaceRefTensorSpace.dim(),-1))

		# Function for transferring fluid stress to structure. Tensor space in reference.
		# self.fluidInterfaceStress = Function(self.interfaceFluidTensorSpace)
		# self.fluidInterfaceStress.vector().zero()

		# traction computation, needs vector space, test and trial functions
		self.tractionFluidVectorSpace = VectorFunctionSpace(meshSolver.mesh, "CG", 1)
		self.tractionTrial = TrialFunction(self.tractionFluidVectorSpace)
		self.tractionTest = TestFunction(self.tractionFluidVectorSpace)

		self.sigma_FSI = Function(self.tensorSpaceRef)

		# to compute stress on reference domain, have dummy velocity and pressure,
		# previous and current time steps.
		self.V = VectorFunctionSpace(meshSolver.mesh, "CG", 2)
		self.Q = FunctionSpace(meshSolver.mesh, "CG", 1)
		self.U_F0 = Function(self.V)
		self.U_F1 = Function(self.V)
		self.P_F0 = Function(self.Q)
		self.P_F1 = Function(self.Q)

		self.ref_ds = Measure('ds', domain = p_s.meshRef, subdomain_data = p_s.facetsRef)
		self.d_FSI = self.ref_ds(20)
		self.n_ref  = FacetNormal(p_s.meshRef)

		# Trial and test functions
		# Define trial and test functions
		self.u = TrialFunction(self.vectorSpace)
		self.p = TrialFunction(self.scalarSpace)
		self.v = TestFunction(self.vectorSpace)
		self.q = TestFunction(self.scalarSpace)

		# Define coefficients
		self.k = Constant(p_s.dt)
		self.U  = 0.5*(self.u0 + self.u)
		self.n  = FacetNormal(self.mesh)
		p_s.mu_f = Constant(p_s.mu_f)
		p_s.rho_f = Constant(p_s.rho_f)


		self.cells = CellFunction("size_t", self.mesh)
		self.facets = FacetFunction("size_t", self.mesh)
		self.dx = Measure('dx', domain = self.mesh, subdomain_data = self.cells)
		self.ds = Measure('ds', domain = self.mesh, subdomain_data = self.facets)

		## beta is a flag. set to zero when periodic boundary conditions are used.
		#self.beta = 1

		# Define symmetric gradient
		def epsilon(u):
		    return sym(nabla_grad(u))

		# Define stress tensor
		def sigma(u, p):
		    return 2*p_s.mu_f*epsilon(u) - p*Identity(len(u))

		# mesh velocity is on degree 1 function space in undeformed domain.
		# copy degree 1 to degree 1
		self.meshVelocityCurrent1.vector()[:] = meshSolver.meshVelocity.vector().get_local()
		# Interpolate or project degree 1 to degree 2.
		self.meshLocal = project(self.meshVelocityCurrent1, self.vectorSpace, solver_type = "mumps", \
			form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )


		# self.sigma_FSI = project(self.sigma_FSI1, self.tensorSpaceRef, solver_type = "mumps", \
			# form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )


		# Tentative velocity step -where to place -self.meshVelocity??
		self.F1 = p_s.rho_f*dot((self.u - self.u0) / self.k, self.v)*self.dx \
				+ p_s.rho_f*dot(dot((self.u0-self.meshLocal), nabla_grad(self.u0)), self.v)*self.dx \
				+ inner(sigma(self.U, self.p0), epsilon(self.v))*self.dx \
				+ dot(self.p0*self.n, self.v)*self.ds \
				- dot(p_s.mu_f*nabla_grad(self.U)*self.n, self.v)*self.ds \
				- dot(self.f, self.v)*self.dx
		self.a1 = lhs(self.F1)
		self.L1 = rhs(self.F1)

		# Pressure correction
		self.a2 = dot(nabla_grad(self.p), nabla_grad(self.q))*self.dx
		self.L2 = dot(nabla_grad(self.p0), nabla_grad(self.q))*self.dx \
				- (1.0/self.k)*div(self.u1)*self.q*self.dx
			#-div(self.u1)*self.q*self.dx
			#-(1.0/self.k)*div(self.u1)*self.q*self.dx

		# Velocity correction
		self.a3 = dot(self.u, self.v)*self.dx
		self.L3 = dot(self.u1, self.v)*self.dx\
				- self.k*dot(nabla_grad(self.p1-self.p0), self.v)*self.dx

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

	def solveFluidProblem(self, p_s, meshInterfaceVelocity, meshDisplacement):
		if self.solverMethod == "IPCS":
			self.IPCS_Fluid_Solver(p_s, meshInterfaceVelocity, meshDisplacement)

		else:
			print "Error. The solver for the fluid problem should be set to IPCS"

	def IPCS_Fluid_Solver(self, p_s, meshInterfaceVelocity, meshDisplacement):
		print ""
		print ""
		print "ENTERING FLUID SOLVER"
		print ''
		print ''

		self.fluidInterfaceVelocity.vector()[:] = meshInterfaceVelocity.vector().get_local()

# meshSolver.meshInterfaceVelocity.vector().get_local()
		# This should occure in initialization....  but requires defined boundaries ...

		[bc.apply(self.A1) for bc in self.bcu]
		[bc.apply(self.A2) for bc in self.bcp]
		[bc.apply(self.A3) for bc in self.bcu]

		### solver.
		I = Identity(2) # could make more generic.

		# Compute tentative velocity step
		#begin("Computing tentative velocity")
		self.b1 = assemble(self.L1)
		[bc.apply(self.b1) for bc in self.bcu]
		solve(self.A1, self.u1.vector(), self.b1, 'bicgstab', 'hypre_amg')
		end()

		# Pressure correction
		#begin("Computing pressure correction")
		self.b2 = assemble(self.L2)
		[bc.apply(self.b2) for bc in self.bcp]
		solve(self.A2, self.p1.vector(), self.b2, 'bicgstab', 'hypre_amg')
		end()

		# Velocity correction
		#begin("Computing velocity correction")
		self.b3 = assemble(self.L3)
		[bc.apply(self.b3) for bc in self.bcu]
		solve(self.A3, self.u1.vector(), self.b3, 'cg', 'sor')

		# Step fluid solution back to reference domain.
		# Calculate stress, then traction
		# Step fluid solution back to reference domain
		self.U_F0.vector()[:] = self.u0.vector()[:]
		self.U_F1.vector()[:] = self.u1.vector()[:]
		self.P_F0.vector()[:] = self.p0.vector()[:]
		self.P_F1.vector()[:] = self.p1.vector()[:]
		# take midpoint of timestep
		self.U_F = 0.5 * (self.U_F0 + self.U_F1)
		self.P_F = 0.5 * (self.P_F0 + self.P_F1)
		# compute stress
		self.tau = p_s.mu_f*(grad(self.U_F) + grad(self.U_F).T)
		self.sigma_FSI1 = -self.P_F*I + self.tau
		# perform piola transform
		self.F = meshDisplacement
		self.DF = grad(self.F) + Identity(2)
		self.J = det(self.DF)

		#self.tau = p_s.mu_f*(grad(self.U_F) + grad(self.U_F).T)
		#self.sigma_FSI1 = -self.P_F*I + self.tau

		# try this sigma (from cbc) that accounts for deformation somehow.
		self.sigma_FSI1 = p_s.mu_f*(grad(self.U_F)*inv(self.DF) + inv(self.DF).T \
		* grad(self.U_F).T) - self.P_F*I

		# Piola transform
		self.sigma_FSI = self.J*self.sigma_FSI1*inv(self.DF).T

		#### have changed name from sigma_FSI1 to sigma_FSI.
		# project onto tensor function space in reference domain,
		# is this necessary?
		# send this to structure.
		# should accoun
		#self.sigma_FSI = project(self.sigma_FSI, self.tensorSpaceRef, solver_type = "mumps", \
			#form_compiler_parameters = {"cpp_optimize" : True, "representation" : "uflacs"} )

		a_F = dot(self.tractionTest, self.tractionTrial)*self.d_FSI
		L_F = -dot(self.tractionTest, dot(self.sigma_FSI, self.n_ref))*self.d_FSI

		# traction is negated here to differ be equal and opposite fluid traction
		A_F = assemble(a_F, keep_diagonal = True)
		A_F.ident_zeros()
		# # , exterior_facet_domains=fluidSolver.fsi_boundary_F)
		B_F = assemble(L_F)

		self.traction = Function(self.tractionFluidVectorSpace)
		solve(A_F, self.traction.vector(), B_F)


		### Debug fluid traction
		# slight discrepency here.
		form = dot(dot(self.sigma_FSI, self.n_ref),self.n_ref)*self.d_FSI
		self.integral_0 = assemble(form)

		# Compute integral of projected (and negated) normal traction
		form = dot(self.traction, self.n_ref)*self.d_FSI
		self.integral_1 = -assemble(form)


		print ""
		print ""
		print "EXITING FLUID SOLVER"
