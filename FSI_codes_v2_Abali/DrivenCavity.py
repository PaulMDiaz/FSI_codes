
from dolfin import *
import numpy as np
import pylab as plt

class DrivenCavity:

	def __init__(self):


                ############## INPUT DATA PARAMETERS ###################

                # Physical parameters
                self.nu_f = 0.01	# Fluid viscosity
                self.nu_s = 0.2	# Structure Poisson coefficient
                self.E_s = 1e3	# Structure Young modulus
                self.rho_f = 1	# Fluid density (incorporated in the fluid corrected pressure as p_corr = p/rho)

                # Numerical parameters
                self.dt = 0.05	# Time step
                self.T = 0.15		#  Set final time for iteration
                self.N = 32		# Number of discretizations (square mesh)

                # Geometric parameters
                self.h = 0.5	# Nondimensional structure thickness
                # Check if N is a multiple of 1/h -> Error check to be included

		# Lame' constants
		self.mu_s = self.E_s/(2.0*(1.0 + self.nu_s))
		self.lambda_s = self.E_s*self.nu_s/((1.0 + self.nu_s)*(1.0 - 2.0*self.nu_s))

		self.mu_f = self.rho_f*self.nu_f

		# Set up a variable for time
		self.t = 0
		self.iter_t = 1

		# Set up the top velocity boundary condition
		class FSS(Expression):

			def eval(self, values, x):

        			if between(x[0], (0.0, 0.3)) and near(x[1], 2.0):
            				values[0] = 0.5*sin(pi*x[0]/0.6)**2
					values[1] = 0.0
        			elif between(x[0], (0.3, 1.7)) and near(x[1], 2.0):
            				values[0] = 0.5
					values[1] = 0.0
        			elif between(x[0], (1.7, 2.0)) and near(x[1], 2.0):
            				values[0] = 0.5*sin(pi*(x[0]-1)/0.6)**2
					values[1] = 0.0
				else:
	    				values[0] = 0.0
					values[1] = 0.0

			def value_shape(self):
				return (2,)

		self.U = FSS(degree = 2)	# Top velocity

		################ DEFINE MESHES AND DOMAINS #######################

		self.mesh = RectangleMesh(Point(0.0, 0.0), Point(2.0, 2.0), self.N, self.N)	# Global mesh
		self.Define_Subdomains()		# Sets the subdomains and the submeshes for fluid and structure

		self.Dim = self.mesh.topology().dim()

		# Variables to generate files
		pwd = './Results_Driven_Cavity_FSI/'
		self.file_d_s = File(pwd + 'd_s.pvd')
		self.file_v_f = File(pwd + 'v_f.pvd')
		self.file_p_f = File(pwd + 'p_f.pvd')

	def Define_Subdomains(self):
		h = self.h
		# Define fluid subdomain of cavity
		class Fluid(SubDomain):
			# Fluid domain is 0 < x < 2.0 and h < y < 2
			def inside(self, x, on_boundary):
				return (between(x[0], (0.0, 2.0)) and between(x[1], (h , 2.0)))
				#return (between(x[0], (0.0, 2.0)) and between(x[1], (h , 2.0)))
		# Define structure subdomain of cavity
		class Structure(SubDomain):
			# Structure domain is 0 < x < 2.0 and 0 < y < h
			def inside(self, x, on_boundary):
				return (between(x[0], (0.0, 2.0)) and between(x[1], (0.0, h)))


		# Initialize interior of entire domain
		# cell function used to mark domains defined by classes above
		self.subdomains = CellFunction('size_t', self.mesh)
		self.subdomains.set_all(0) 	# Set entire domain to 0

		# Initialize classes for fluid and structure domains
		fluid = Fluid()
		structure = Structure()

		# Mark fluid and structure domains
		# fluid.mark(interior, 0)    Already marked as 0
		structure.mark(self.subdomains, 1)

		# Define submeshes for fluid and structure on subdomains
		self.mesh_f = SubMesh(self.mesh, self.subdomains, 0)
		self.mesh_s = SubMesh(self.mesh, self.subdomains, 1)

	def Define_Boundary_Conditions(self, S, F):
		# def Define_Boundaries(self, Structure, Fluid, FSI_Boundary):
		#	Description: Defines domains for fluid and structure solvers
		#				 and sets boundary conditions.
		#   Input: Structure - structure solver object
		#		   Fluid - fluid solver object
		#		   FSI_Boundary - thickness of bottom wall


		############## Define boundary domain locations #################
		h = self.h

		# Top boundary for free-stream fluid flow
		class Top(SubDomain):
			def inside(self, x, on_boundary):
				return near(x[1], 1.0)

		# Left boundary for structure and fluid
		class Left(SubDomain):
			def inside(self, x, on_boundary):
				return near(x[0], 0.0)

		# Right boundary for structure and fluid
		class Right(SubDomain):
			def inside(self, x, on_boundary):
				return near(x[0], 1.0)

		# FSI boundary for fluid and structure
		class FSI(SubDomain):
			def inside(self, x, on_boundary):
				return near(x[1], h)

		# Bottom boundary for structure
		class Bottom(SubDomain):
			def inside(self, x, on_boundary):
				return near(x[1], 0.0)

		############# Initialize structure boundary #############

		# Create mesh function for Structure boundary and mark all 0
		S.cells = CellFunction("size_t", S.mesh)
		S.facets = FacetFunction("size_t", S.mesh)
		S.facets.set_all(0)

		# Initialize boundary objects
		S.left = Left()
		S.right = Right()
		S.fsi = FSI()
		S.bottom = Bottom()

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

		# Create fluid boundary objects
		F.left = Left()
		F.right = Right()
		F.fsi = FSI()
		F.top = Top()

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
		noSlipLeft = DirichletBC(F.V_space, Constant((0, 0)), F.left)
		noSlipRight = DirichletBC(F.V_space, Constant((0, 0)), F.right)
		#  Freestream velocity boundary condition for top of cavity
		freestreamV = DirichletBC(F.V_space, self.U, F.top)
		# Initialized as zero, equal to the initial velocity field
		fluidFSI = DirichletBC(F.V_space, F.u_mesh, F.fsi)
		# Pressure
		F.bcp = [DirichletBC(F.S_space, Constant(0.0), F.left), \
					DirichletBC(F.S_space, Constant(0.0), F.right), \
					DirichletBC(F.S_space, Constant(0.0), F.top)]

		# Set up the boundary conditions
		F.bcu = [noSlipLeft, noSlipRight, freestreamV, fluidFSI]

	def Save_Results(self, S, F):

		S.d_.assign(S.d)
		self.file_d_s << (S.d_, self.t)
		F.v_.assign(F.u1)
		self.file_v_f << (F.v_, self.t)
		F.p_.assign(F.p1)
		self.file_p_f << (F.p_, self.t)
