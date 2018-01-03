
from dolfin import *
import numpy as np
import pylab as plt
from mshr import *

# Construct a structure solver on its own. Solve for both displacement and Velocity
#

set_log_level(ERROR)
ElementType = "CG" # CG is equivilent to lagrange.

#  Set the Structure Element Degree
ElementDegree = 1 # Standard linear lagrange element.

#B = Constant((0.0, -1.0))
B = Constant((0.0, 0.0))

########################

#nu_f = 0.001	# Fluid viscosity (was 0.2)
nu_s = 0.4	# Structure Poisson coefficient should be 0.2
mu_s = 0.5e6 # structure first lame constant (very small)

#rho_f = 1000.0	# Fluid density (incorporated in the fluid corrected pressure as p_corr = p/rho)

rho_s = 1000.0

E_s = mu_s*2*(1+nu_s)	# structure elastic modulus
lambda_s = nu_s*E_s/((1+nu_s)*(1-2*nu_s)) # 2nd lame constant

# mean inlet velocity
#U_mean = 0.2

# Numerical parameters
dt = 0.001	# Time step
#T = 10.00	#  Set final time for iteration
T = 3.000	#  Set final time for iteration
N = 64		# Number of discretizations (square mesh) (place cell edge on FSI)

# Geometric parameters
L = 2.5 	#Channel length
H = 0.41	#Channel height
l = 0.35 	#Bar length
h = 0.02	#Bar height

#s_b = 2*h 	#Border about structure over which to refine.

# x coordinate of start of bar
x_bar = 0.6-l

#mu_f = rho_f*nu_f

# Set up a variable for time
t_vec = np.linspace(0.0,T,T/dt+1)
t = t_vec[0]
iter_t = 1

################ DEFINE MESHES AND DOMAINS #######################

# this definition leaves bar tangent to circle with small gap
bar = Rectangle(Point(x_bar,0.19), Point(0.6, 0.19+h))
# Ensures overlap between bar and circle.
#bar_2 = Rectangle(Point(x_bar-0.01,0.19), Point(0.6, 0.19+h))

# Monitor point for later use
A_point = Point(0.6, 0.2)

# structure domain
s_domain = bar

# total domain

# generate global mesh
mesh = generate_mesh(s_domain, N)



# Variational spaces
S_space = FunctionSpace(mesh, ElementType, ElementDegree)
V_space = VectorFunctionSpace(mesh, ElementType, ElementDegree)
T_space = TensorFunctionSpace(mesh, ElementType, ElementDegree)

# Mixed space for computing both displacement and rate of change of displacement
V_element = VectorElement(ElementType, mesh.ufl_cell() , ElementDegree)
#V2_temp = V_element*V_element

#V2_space = FunctionSpace(mesh, V2_temp)
V2_space = FunctionSpace(mesh, MixedElement([V_element, V_element]))
V2_space = FunctionSpace(mesh, V_element*V_element)

V = TestFunction(V2_space)
dU = TrialFunction(V2_space)
U = Function(V2_space)
U0 = Function(V2_space)

# Load initial conditions to u0 and v0. Otherwise set to 0.
u0 = Constant((0,)*V_space.mesh().geometry().dim())
v0 = Constant((0,)*V_space.mesh().geometry().dim())

# Functions for solver
u_t, v_t = split(V) 	# Test functions
u, v = split(U)		# Functions

cells = CellFunction("size_t", mesh)
dx = Measure('dx', domain = mesh, subdomain_data = cells)

# Project u0 and v0 into U0
a_proj = inner(dU, V)*dx
L_proj = inner(u0, u_t)*dx + inner(v0, v_t)*dx
solve(a_proj == L_proj, U0)

# load data
#U_temp = np.loadtxt('nodal_U_64')
#U0.vector()[:] = U_temp

u0, v0 = split(U0)


## Function for results
u_res = Function(V_space, name = 'u')
v_res = Function(V_space, name = 'v')



Dim = mesh.topology().dim()

# Variables to generate files
pwd = './Results_Cylinder_bar_T/'
file_u_s = File(pwd + 'u_s.pvd')
file_v_s = File(pwd + 'v_s.pvd')
#file_v_f = File(pwd + 'v_f.pvd')
#file_p_f = File(pwd + 'p_f.pvd')

def compute_forces(self,Mesh,nu,u,p, ds):
	mesh = Mesh
	nu = nu
	# should include both pressure contribution and shear forces.
	# Face normals
	n = FacetNormal(mesh)
	# Stress tensor
	# Traction
	sigma = nu*(grad(u)+grad(u).T) -p*Identity(2)
	T = dot(sigma,n)

	drag = -T[0]*ds(3)-T[0]*ds(4)
	lift = -T[1]*ds(3)-T[1]*ds(4)

	drag1 = assemble(drag); lift1 = assemble(lift);
	return drag1, lift1

############## Define boundary domain locations #################

class Fsi(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > x_bar - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > DOLFIN_EPS and x[1] < (0.19 + DOLFIN_EPS) and on_boundary or \
		x[0] > x_bar - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > (0.21 - DOLFIN_EPS) and x[1] < H - DOLFIN_EPS and on_boundary or \
		x[0] > 0.6 - DOLFIN_EPS and x[0] < L - DOLFIN_EPS and \
        x[1] > DOLFIN_EPS and x[1] < H - DOLFIN_EPS and on_boundary

#cylinder = CompiledSubDomain('on_boundary && x[0] > DOLFIN_EPS && x[0]< x_bar + DOLFIN_EPS && x[1]>DOLFIN_EPS && x[1] < H-DOLFIN_EPS', L = L, H = H, x_bar = x_bar)
#fsi = CompiledSubDomain('on_boundary && x[0] > x_bar + DOLFIN_EPS && x[0]< L - DOLFIN_EPS && x[1]>DOLFIN_EPS && x[1] < H-DOLFIN_EPS', L = L, H = H, x_bar = x_bar)

class Left(SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > x_bar - DOLFIN_EPS and x[0] < x_bar + DOLFIN_EPS and \
        x[1] > DOLFIN_EPS and x[1] < H - DOLFIN_EPS and on_boundary

# Compile subdomains
##walls = Walls()
#inlet = Inlet()
#outlet = Outlet()
#cylinder = Cylinder()
fsi = Fsi()
left = Left()

#inlet   = CompiledSubDomain('near(x[0], 0) && on_boundary ')
#outlet  = CompiledSubDomain('near(x[0], L) && on_boundary', L = L)
#walls   = CompiledSubDomain('near(x[1], 0) || near(x[1], H) && on_boundary', H = H)

############# Initialize structure boundary #############

# Create mesh function for Structure boundary and mark all 0
cells = CellFunction("size_t", mesh)
facets = FacetFunction("size_t", mesh)
facets.set_all(0)

#bottom = CompiledSubDomain('x[1] < DOLFIN_EPS && on_boundary')

############# Initialize structure boundary #############
# Initialize boundary objects

# Mark boundaries
left.mark(facets, 1)
fsi.mark(facets, 2)

ds = Measure('ds', domain = mesh, subdomain_data = facets)
dx = Measure('dx', domain = mesh, subdomain_data = cells)

n = FacetNormal(mesh)

#  BCs for the left side (no displacement)

LeftBC = DirichletBC(V2_space, Constant((0.0, 0.0, 0.0, 0.0)), left)

#  Set up the boundary conditions
bcs = [LeftBC]

############ Initialize fluid boundary ###################


count = 0

while t <= T: # + DOLFIN_EPS:
	print 'STARTING TIME LOOP ITERATION ', count

	t += dt

	for ii in range(3): #change to 3 to run properly.

		# displacements and velocities at mid points
		u_mid = 0.5*(u0 + u)
		v_mid = 0.5*(v0 + v)

		I = Identity(Dim)

		## Project stress from fluid tensor space onto structure tensor space
		#sigma_FSI = project(F.sigma_FSI, T_space, solver_type = "mumps",\
		#form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

		# Kinematics
		# Only use 2nd pk stres evaluated at u_mid it seems.
		F =  variable(I + grad(u_mid))		# Deformation gradient tensor
		C =  variable(F.T*F)			# Right Cauchy-Green tensor

		# Green lagrange stress tensor (green-st venant strain tensor)
		E = variable(0.5*(C-I))

		# stored strain energy
		psi = lambda_s/2*(tr(E)**2)+mu_s*tr(E*E)

		# 2nd piola kirchoff stress tensor
		S2 = diff(psi, E)

		# 1st piola kirchoff stress tensor
		#evaluated at u_mid
		S1 = F*S2

		# Invariants of deformation tensor
		J = det(F)

		# Evaluate

		k = Constant(dt)
		# Variational form for hyperelasticity

		# Traction on boundary

		# solves and gives 0 displacment
		#T_hat = Constant((0.0, 0.0))

		# Solves and gives too great a displacement.
		n_s = FacetNormal(mesh)

		# To be consistent, use u_mid for T_hat too...
		#T_hat = (J*inv(F)*sigma_FSI).T * n_s
		# gives same results as index method
		T_hat = Constant((1.0, -10.0))

		# The variational form corresponding to hyperelasticity

		L = rho_s*inner(v - v0, u_t)*dx \
		+ k*inner(S1, grad(u_t))*dx \
		- k*inner(B, u_t)*dx \
		+ inner(u - u0, v_t)*dx \
		- k*inner(v_mid, v_t)*dx

		L = L-k*inner(T_hat, u_t)*ds(2)


		# introduced in line 499 of solution_algorithms
		a = derivative(L, U, dU)

		# Setup problem
		# may have to fix bcu to align with Twist.
		# a bit iffy on how solve, step, update interact in CBC twist... write out, comb over.

		problem = NonlinearVariationalProblem(L, U, bcs, a)
		solver = NonlinearVariationalSolver(problem)
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["relative_tolerance"] = 1e-12
		solver.parameters["newton_solver"]["maximum_iterations"] = 100

		solver.solve()

	# propogate displacements and velocities
	U0.assign(U)

	u2, v2 = U.split()
	# maybe like this, or maybe save U...
	file_u_s << (u2, t)
	file_v_s << (v2, t)

	count += 1
