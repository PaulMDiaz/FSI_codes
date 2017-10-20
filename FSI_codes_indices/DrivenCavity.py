
from dolfin import *
import numpy as np
import pylab as plt
from mshr import *

# edit to make similar to chorin that works. Remove pressure BC. Rearrange fluid solver to match. Comment out previous arrangement.
# breaks after a couple of iterations ' unable to solve linear system using PETSC Krylov solver' solution failed to converge.
# breaks when mesh nodes do not align perfectly with interface. If number of discretizations is selected carefully then the code runs. Evidently
#this is a working and not elegant or long term solution. Possibly if the mesh deforms the same error will reappear.
# Therefore mesh still has to be fixed for long term solution such that fluid and structure boundaries align perfectly with mesh nodes.

########################

# when I run it with linear and E_s = 1e3. It breaks at loop iteration time 0.2.
# Mesh set to the generate_mesh  method. Norms are 8 e-08 and 4e-18.

# In the other method: Still breaks at loop iteration time 0.2. first norm is 3e-9. Marginally different.

# Ignoring precision in integral metadata compiled using quadrature representation. Not implemented.

class DrivenCavity:

	def __init__(self):
        ############## INPUT DATA PARAMETERS ###################
	    # Physical parameters
		self.nu_f = 0.2	# Fluid viscosity (was 0.2)
		self.nu_s = 0.2	# Structure Poisson coefficient should be 0.2
		self.E_s = 1e3	# Structure Young modulus (was 1e3)
		self.rho_f = 1.0	# Fluid density (incorporated in the fluid corrected pressure as p_corr = p/rho)
		self.rho_s = 1.0

		#self.rho_s = 1

		# Numerical parameters
		self.dt = 0.01	# Time step
		self.T = 0.03	#  Set final time for iteration
		self.N = 64		# Number of discretizations (square mesh) (place cell edge on FSI)

		# Geometric parameters
		self.h = 0.5	# Nondimensional structure thickness
		self.H = 2.0	# Height of entire domain (1.5 for validation, should be 2)
		self.W = 2.0		# Width of entire domain (1 for validation, should be 2)
		# Check if N is a multiple of 1/h -> Error check to be included

		# Lame' constants
		# first lame constant is the shear modulus
		# should be approx 2e6
		self.mu_s = Constant(self.E_s/(2.0*(1.0 + self.nu_s)))

		# second lame constant is bulk modulus less 2/3 shear modulus
		self.lambda_s = self.E_s*self.nu_s/((1.0 + self.nu_s)*(1.0 - 2.0*self.nu_s))

		self.mu_f = self.rho_f*self.nu_f

		# Set up a variable for time
		self.t = 0
		self.iter_t = 1

		################ DEFINE MESHES AND DOMAINS #######################

		# This fixes the issue of nodes falling on the interface.
		# It makes the subsequent subdomain definitions appear cumbersome.
		# Perhaps these can be eliminated?
		# Worried that these domains won't match with mesh and fluid domains once deformation happens...


		#self.domain = Rectangle(Point(0.0, 0.0), Point(self.W, self.H))
		#self.f_domain = Rectangle(Point(0.0, self.h), Point(self.W, self.H))
		#self.s_domain = Rectangle(Point(0.0, 0.0), Point(self.W, self.h))

		#self.domain.set_subdomain(1,self.s_domain)
		#Interface = Line()

		self.mesh = RectangleMesh(Point(0.0, 0.0), Point(self.W, self.H), self.N, self.N)	# Global mesh

		#self.mesh = generate_mesh(self.domain, self.N)
		# maybe use snap_boundary()
		self.Define_Subdomains()		# Sets the subdomains and the submeshes for fluid and structure

		self.Dim = self.mesh.topology().dim()

		# Variables to generate files
		pwd = './Results_Driven_Cavity_FSI/'
		self.file_d_s = File(pwd + 'd_s.pvd')
		self.file_v_f = File(pwd + 'v_f.pvd')
		self.file_p_f = File(pwd + 'p_f.pvd')

	def Define_Subdomains(self):
		h = self.h
		W = self.W
		H = self.H

		# Fails to converge: ##
		#Define tolerance for subdomains. (as in tutorial pg 97)
		tol = 1E-14

		# Define fluid subdomain of cavity
		class Fluid(SubDomain):
			### Fluid domain is 0 < x < 2.0 and h < y < 2
			def inside(self, x, on_boundary):
				###return True if 0.0 <= x[0] <= 2.0 and  h <= x[1] <=  2.0 else False
				##return (between(x[0], (0.0, W)) and between(x[1], (h , H)))
				return True if x[1] >= h - tol else False
		### Define structure subdomain of cavity
		class Structure(SubDomain):
			### Structure domain is 0 < x < 2.0 and 0 < y < h
			def inside(self, x, on_boundary):
				#return True if 0.0 <= x[0] <= 2.0 and 0.0 <= x[1] <=  h else False
				##return (between(x[0], (0.0, W)) and between(x[1], (0.0, h)))
				return True if x[1] <= h + tol else False

		## Replace with
		#structure = CompiledSubDomain('x[0] >= 0.0 && x[0] <= W1 && x[1] <= h1', h1 = h, W1 = W)
		#fluid = CompiledSubDomain('x[0] >= 0.0 && x[0] <= W1 && x[1] >= h1 && x[1] <= H1', W1 = W, h1 = h, H1 = H)

		### Initialize interior of entire domain
		### cell function used to mark domains defined by classes above
		self.subdomains = CellFunction('size_t', self.mesh)
		# was CellFunction
		self.subdomains.set_all(0) 	# Set entire domain to 0

		## Initialize classes for fluid and structure domains
		fluid = Fluid()
		structure = Structure()


		# Mark fluid and structure domains
		fluid.mark(self.subdomains, 0)
		structure.mark(self.subdomains, 1)

		# Define submeshes for fluid and structure on subdomains
		# This might be a mistake. Doesn't happen in fenics tutorial pg 100. Try without.
		self.mesh_f = SubMesh(self.mesh, self.subdomains, 0)
		self.mesh_s = SubMesh(self.mesh, self.subdomains, 1)


		# Set up the top velocity boundary condition
		class FSS(Expression):
			def eval(self, values, x):
				if between(x[0], (0.0, 0.3)) and near(x[1], 2):
					values[0] = 0.5*sin(pi*x[0]/0.6)**2
					values[1] = 0.0
				elif between(x[0], (0.3, 1.7)) and near(x[1], 2):
					values[0] = 0.5
					values[1] = 0.0
				else:
					values[0] = 0.5*sin(pi*(x[0]-2)/0.6)**2
					values[1] = 0.0
					# elif between(x[0], (1.7, H)) and near(x[1], H)::

			# value shape is necessary to make expression vector or tensor valued
			def value_shape(self):
				return (2,)

		self.V = VectorFunctionSpace(self.mesh_f, 'P', 2)

		#self.U_top = FSS(degree = 1)	# Top velocity
		self.U_top = interpolate(FSS(degree = 0), self.V)
		#self.U_top = Expression(("1" , "0.0"), degree = 0)	# Top velocity

	def Define_Boundary_Conditions(self, S, F):
		# def Define_Boundaries(self, Structure, Fluid, FSI_Boundary):
		#	Description: Defines domains for fluid and structure solvers
		#				 and sets boundary conditions.
		#   Input: Structure - structure solver object
		#		   Fluid - fluid solver object
		#		   FSI_Boundary - thickness of bottom wall


		############## Define boundary domain locations #################
		h = self.h
		W = self.W
		H = self.H

		# Top boundary for free-stream fluid flow
		#class Top(SubDomain):
			#def inside(self, x, on_boundary):
				#return near(x[1], H)

		## Left boundary for structure and fluid
		#class Left(SubDomain):
			#def inside(self, x, on_boundary):
				#return near(x[0], 0.0)

		# Right boundary for structure and fluid
		#class Right(SubDomain):
			#def inside(self, x, on_boundary):
				#return near(x[0], W)

		# FSI boundary for fluid and structure
		#class FSI(SubDomain):
			#def inside(self, x, on_boundary):
				#return near(x[1], h)

		 #Bottom boundary for structure
		#class Bottom(SubDomain):
			#def inside(self, x, on_boundary):
				#return near(x[1], 0.0)

		############# Initialize structure boundary #############

		# Create mesh function for Structure boundary and mark all 0
		S.cells = CellFunction("size_t", S.mesh)
		S.facets = FacetFunction("size_t", S.mesh)
		S.facets.set_all(0)

		top = CompiledSubDomain('near(x[1], H1) && on_boundary', H1 = H)
		#top = CompiledSubDomain('x[1] > H1 - DOLFIN_EPS && on_boundary', H1 = H)
		# Left boundary for structure and fluid
		left = CompiledSubDomain('near(x[0], 0.0) && on_boundary')
		#left = CompiledSubDomain('x[0] < DOLFIN_EPS && on_boundary')
		# Right boundary for structure and fluid
		#right = CompiledSubDomain('x[0] > W1 - DOLFIN_EPS && on_boundary', W1 = W)
		right = CompiledSubDomain('near(x[0],W1) && on_boundary', W1 = W)

		# FSI boundary for fluid and structure
		#fsi = CompiledSubDomain('x[1] > h1 - DOLFIN_EPS || x[1] < h1+ DOLFIN_EPS && on_boundary', h1 = h)
		fsi = CompiledSubDomain('near(x[1], h1)  && on_boundary', h1 = h)
		# Bottom boundary for structure
		bottom = CompiledSubDomain('near(x[1], 0.0) && on_boundary')
		#bottom = CompiledSubDomain('x[1] < DOLFIN_EPS && on_boundary')

		############# Initialize structure boundary #############
		# Initialize boundary objects
		S.left = left
		S.right = right
		S.fsi = fsi
		S.bottom = bottom

		## Initialize boundary objects
		#S.left = Left()
		#S.right = Right()
		#S.fsi = FSI()
		#S.bottom = Bottom()

		# Mark boundaries
		S.left.mark(S.facets, 1)
		S.right.mark(S.facets, 2)
		S.fsi.mark(S.facets, 3)
		# S.bottom.mark(S.facets, 0)

		S.dA = Measure('ds', domain = S.mesh, subdomain_data = S.facets)
		S.dV = Measure('dx', domain = S.mesh, subdomain_data = S.cells)

		S.N = FacetNormal(S.mesh)

		#  BCs for the left side (no displacement)
		LeftBC = DirichletBC(S.V_space, Constant((0.0, 0.0)), S.left)

		#  BCs for the right side (no displacement)
		RightBC = DirichletBC(S.V_space, Constant((0.0, 0.0)), S.right)

		# BC for the bottom wall (no displacement) # it can be changed to no displacement variation by deleting this B.C.
		BottomBC = DirichletBC(S.V_space, Constant((0.0, 0.0)), S.bottom)

		#  Set up the boundary conditions
		S.bcs = [LeftBC, RightBC, BottomBC]

		############ Initialize fluid boundary ###################

		# Create mesh function and mark all 0
		F.cells = CellFunction("size_t", F.mesh)
		F.facets = FacetFunction("size_t", F.mesh)
		F.facets.set_all(0)

		#interaction.mark(S.facets,2)

		############ Initialize fluid boundary ###################

		# Initialize boundary objects
		F.left = left
		F.right = right
		F.fsi = fsi
		F.top = top

		# Create fluid boundary objects
		#F.left = Left()
		#F.right = Right()
		#F.fsi = FSI()
		#F.top = Top()

		# mark boundaries for fluid
		F.left.mark(F.facets, 1)
		F.right.mark(F.facets, 2)
		F.fsi.mark(F.facets, 3)
		# F.top.mark(F.facets, 0)

		F.da = Measure('ds', domain = F.mesh, subdomain_data = F.facets)
		F.dx = Measure('dx', domain = F.mesh, subdomain_data = F.cells)

		F.n = FacetNormal(F.mesh)

		#  Define  fluid boundary conditions
		#  Noslip boundary condition for bottom and walls of cavity
		noSlipLeft = DirichletBC(F.V_space, Constant((0.0, 0.0)), F.left)
		noSlipRight = DirichletBC(F.V_space, Constant((0.0, 0.0)), F.right)
		#  Freestream velocity boundary condition for top of cavity
		# run FSI problem:
		freestreamV = DirichletBC(F.V_space, self.U_top, F.top)
		#run comparison to ghia et al
		#freestreamV = DirichletBC(F.V_space, Constant((1.0,0)), F.top)
		# Initialized as zero, equal to the initial velocity field
		#fluidFSI = DirichletBC(F.V_space, Constant((0, 0)), F.fsi)
		fluidFSI = DirichletBC(F.V_space, F.u_mesh, F.fsi)

		# Pressure
		#F.bcp = [DirichletBC(F.S_space, Constant(0.0), F.left), DirichletBC(F.S_space, Constant(0.0), F.right), DirichletBC(F.S_space, Constant(0.0), F.top)]

		# Set up the boundary conditions
		F.bcu = [noSlipLeft, noSlipRight, freestreamV, fluidFSI]

	def Save_Results(self, S, F):

		S.d_.assign(S.d)
		self.file_d_s << (S.d_, self.t)
		F.v_.assign(F.u1)
		self.file_v_f << (F.v_, self.t)
		F.p_.assign(F.p1)
		self.file_p_f << (F.p_, self.t)
