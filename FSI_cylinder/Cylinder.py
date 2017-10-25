
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

# Problem specific section. ie lid driven cavity, or cylinder... have class generic?
class problem_specific:

	def __init__(self):
        ############## INPUT DATA PARAMETERS ###################
	    # Physical parameters
		self.nu_f = 0.001	# Fluid viscosity (was 0.2)
		self.nu_s = 0.4	# Structure Poisson coefficient should be 0.2
		self.mu_s = 0.5e6 # structure first lame constant
		#self.E_s = 1e5	# Structure Young modulus (was 1e3)
		self.rho_f = 1000.0	# Fluid density (incorporated in the fluid corrected pressure as p_corr = p/rho)
		self.rho_s = 1000.0

		self.E_s = self.mu_s*2*(1+self.nu_s)	# structure elastic modulus
		self.lambda_s = self.nu_s*self.E_s/((1+self.nu_s)*(1-2*self.nu_s)) # 2nd lame constant

		# mean inlet velocity
		self.U_mean = 0.2

		# Numerical parameters
		self.dt = 0.001	# Time step
		self.T = .005	#  Set final time for iteration
		self.N = 64		# Number of discretizations (square mesh) (place cell edge on FSI)

		# Geometric parameters
		self.L = 2.5 	#Channel length
		self.H = 0.41	#Channel height
		self.l = 0.35 	#Bar length
		self.h = 0.02	#Bar height

		# x coordinate of start of bar
		self.x_bar = 0.6-self.l

		self.mu_f = self.rho_f*self.nu_f

		# Set up a variable for time
		self.t = 0
		self.iter_t = 1

		################ DEFINE MESHES AND DOMAINS #######################

		self.channel = Rectangle(Point(0, 0), Point(self.L, self.H))
		self.cylinder = Circle(Point(0.2, 0.2), 0.05,50)
		# this definition leaves bar tangent to circle with small gap
		self.bar = Rectangle(Point(self.x_bar,0.19), Point(0.6, 0.19+self.h))
		# Ensures overlap between bar and circle.
		self.bar_2 = Rectangle(Point(self.x_bar-0.01,0.19), Point(0.6, 0.19+self.h))

		# Monitor point for later use
		self.A_point = Point(0.6, 0.2)

		# fluid domain
		self.f_domain = self.channel - (self.cylinder + self.bar_2) - (self.bar_2-self.bar)

		# structure domain
		self.s_domain = self.bar

		# total domain
		#self.domain = self.channel - (self.cylinder + self.bar_2) + self.bar_2
		self.domain = self.f_domain +self.s_domain


		# Set structure as subdomain, first type for global mesh
		self.domain.set_subdomain(1,self.s_domain)
		#self.domain.set_subdomain(2,self.f_domain)

		#print self.domain.get_subdomains.__getattribute__
		#print dir(self.domain.get_subdomains.__getattribute__)
		#print dir(self.domain.has_subdomains)
		# generate global mesh
		self.mesh = generate_mesh(self.domain, self.N)


		#self.Define_Subdomains()		# Sets the subdomains and the submeshes for fluid and structure

		#fluid   = CompiledSubDomain(self.f_domain)

		#create a MeshFunction with unsigned integer values (the subdomain numbers)
		# with dimension 2, which is the cell dimension of this problem.
		#markers = MeshFunction('size_t', self.mesh, 2, self.mesh.domains())

		#markers = MeshFunction('size_t', p_s.mesh, 2, p_s.mesh.domains())

		#cell_markers = SubsetIterator(markers, 1)

		#mesh = refine(mesh, cell_markers)

	#def Define_Subdomains(self):
		# Define fluid subdomain
		#f_domain = self.f_domain
		#s_domain = self.s_domain

		#class Fluid(SubDomain):
			#def inside(self, x, on_boundary):
				#return True if in f_domain else False

		# Define structure subdomain
		#class Structure(SubDomain):
			#def inside(self, x, on_boundary):
				#return True if in s_domain else False

		#self.fluid = Fluid()
		#self.structure = Structure()

		# Set fluid and structure as subdomains, type 2 for submesh.
		#self.subdomains = CellFunction('size_t', self.mesh)
		#self.subdomains.set_all(0)

		# Mark fluid and structure domains
		#self.fluid.mark(self.subdomains, 0)
		#self.structure.mark(self.subdomains, 1)

		# Create submeshes
		#self.mesh_f = SubMesh(self.mesh, self.subdomains, 0)
		#self.mesh_s = SubMesh(self.mesh, self.subdomains, 1)

		self.mesh_f = SubMesh(self.mesh, 0)
		self.mesh_s = SubMesh(self.mesh, 1)
		# do I need these subdomains if I already have teh domain.set_subdomain.
		# It looks as though I could do one or the other. Try with just the above.
		#self.subdomains = CellFunction('size_t', self.mesh)
		#self.subdomains.set_all(0)

		## Sets the subdomains and the submeshes for fluid and structure
		#self.Define_Subdomains()

		self.Dim = self.mesh.topology().dim()

		# Variables to generate files
		pwd = './Results_Cylinder_FSI/'
		self.file_d_s = File(pwd + 'd_s.pvd')
		self.file_v_f = File(pwd + 'v_f.pvd')
		self.file_p_f = File(pwd + 'p_f.pvd')


	def Define_Boundary_Conditions(self, S, F):

		L = self.L
		H = self.H
		x_bar = self.x_bar
		U_mean = self.U_mean
		############## Define boundary domain locations #################

		# It would be nice to lump walls and cylinder together.
		inlet   = CompiledSubDomain('near(x[0], 0) && on_boundary ')
		outlet  = CompiledSubDomain('near(x[0], L) && on_boundary', L = L)
		walls   = CompiledSubDomain('near(x[1], 0) || near(x[1], H) && on_boundary', H = H)

		# Distinguish between cylinder support and cantilver beam with fsi surface
		# current formulation of cylinder might obtain not only circle but also weird cutout of rectangle... could be nasty.
		# Should probably edit implementation of Dolfin_eps so that vertices lie on both subdomains..
		cylinder = CompiledSubDomain('on_boundary && x[0] > DOLFIN_EPS && x[0]< x_bar + DOLFIN_EPS && x[1]>DOLFIN_EPS && x[1] < H-DOLFIN_EPS', L = L, H = H, x_bar = x_bar)
		fsi = CompiledSubDomain('on_boundary && x[0] > x_bar + DOLFIN_EPS && x[0]< L - DOLFIN_EPS && x[1]>DOLFIN_EPS && x[1] < H-DOLFIN_EPS', L = L, H = H, x_bar = x_bar)

		# Need to identify LHS of bar.
		# Try within tol of x_bar... Otherwise see Abali
		left = CompiledSubDomain('on_boundary && x[0] > x_bar - DOLFIN_EPS && x[0]< x_bar + DOLFIN_EPS && x[1]>DOLFIN_EPS && x[1] < H-DOLFIN_EPS', L = L, H = H, x_bar = x_bar)
		fsi = CompiledSubDomain('on_boundary && x[0] > x_bar + DOLFIN_EPS && x[0]< L - DOLFIN_EPS && x[1]>DOLFIN_EPS && x[1] < H-DOLFIN_EPS', L = L, H = H, x_bar = x_bar)

		############# Initialize structure boundary #############

		# Create mesh function for Structure boundary and mark all 0
		S.cells = CellFunction("size_t", S.mesh)
		S.facets = FacetFunction("size_t", S.mesh)
		S.facets.set_all(0)

		#bottom = CompiledSubDomain('x[1] < DOLFIN_EPS && on_boundary')

		############# Initialize structure boundary #############
		# Initialize boundary objects
		S.left = left
		S.fsi = fsi

		# Mark boundaries
		S.left.mark(S.facets, 1)
		S.fsi.mark(S.facets, 2)

		S.ds = Measure('ds', domain = S.mesh, subdomain_data = S.facets)
		S.dx = Measure('dx', domain = S.mesh, subdomain_data = S.cells)

		S.n = FacetNormal(S.mesh)

		#  BCs for the left side (no displacement)
		# Structure is in essence a cantilver beam.
		LeftBC = DirichletBC(S.V_space, Constant((0.0, 0.0)), S.left)

		#  Set up the boundary conditions
		S.bcs = [LeftBC]

		# Presently no BC for FSI.

		############ Initialize fluid boundary ###################

		# Create mesh function and mark all 0
		F.cells = CellFunction("size_t", F.mesh)
		F.facets = FacetFunction("size_t", F.mesh)
		F.facets.set_all(0)

		############ Initialize fluid boundary ###################

		# Initialize boundary objects
		F.inlet = inlet
		F.outlet = outlet
		F.fsi = fsi
		F.cylinder = cylinder
		F.walls = walls

		# mark boundaries for fluid
		F.inlet.mark(F.facets, 1)
		F.outlet.mark(F.facets, 2)
		F.fsi.mark(F.facets, 3)
		F.cylinder.mark(F.facets,4)
		F.walls.mark(F.facets,5)

		F.ds = Measure('ds', domain = F.mesh, subdomain_data = F.facets)
		F.dx = Measure('dx', domain = F.mesh, subdomain_data = F.cells)

		F.n = FacetNormal(F.mesh)

		#  Define  fluid boundary conditions
		#  Noslip boundary condition for bottom and walls of cavity
		noSlipwalls = DirichletBC(F.V_space, Constant((0.0, 0.0)), F.walls)
		noSlipcylinder = DirichletBC(F.V_space, Constant((0.0, 0.0)), F.cylinder)
		#  Freestream velocity boundary condition for top of cavity
		# run FSI problem:

		# Define inflow profile
		inlet_profile = ('1.5*U_mean*x[1]*(H-x[1])/pow(H/2, 2)', '0')

		bcu_inlet = DirichletBC(F.V_space, Expression(inlet_profile, H = H, U_mean = U_mean, degree=2), F.inlet)
		bcu_walls = DirichletBC(F.V_space, Constant((0, 0)), F.walls)
		bcu_cylinder = DirichletBC(F.V_space, Constant((0, 0)), F.cylinder)

		bcu_fsi = DirichletBC(F.V_space, Constant((0, 0)), F.fsi)
		#bcu_fsi = DirichletBC(F.V_space, F.u_mesh, F.fsi)

		# Pressure
		bcp_outlet = DirichletBC(F.S_space, Constant(0), F.outlet)

		# Set up the boundary conditions
		F.bcu = [bcu_inlet, bcu_walls, bcu_cylinder, bcu_fsi]
		F.bcp = [bcp_outlet]

# compute and save nodal values
#nodal_values_u = F.u_.vector().array()
#np.savetxt('nodal_u', nodal_values_u)
#nodal_values_p = F.p_.vector().array()
#np.savetxt('nodal_p', nodal_values_p)
#nodal_values_d = S.d_.vector().array()
#np.savetxt('nodal_d', nodal_values_d)

	def Save_Results(self, S, F):

		S.d_res.assign(S.d_)
		self.file_d_s << (S.d_res, self.t)
		F.u_res.assign(F.u_)
		self.file_v_f << (F.u_res, self.t)
		F.p_res.assign(F.p_)
		self.file_p_f << (F.p_res, self.t)
