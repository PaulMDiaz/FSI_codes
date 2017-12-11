"""
FEniCS tutorial demo program: Linear elastic problem.
  -div(sigma(u)) = f
The model is used to simulate an elastic beam clamped at
its left end and deformed under its own weight.
"""

# Put E and nu into this solver...

from __future__ import print_function
from fenics import *
from dolfin import *
import numpy as np
from mshr import *

# Scaled variables
W = 2; h = 0.5; H = 2;
mu = 1                  # lame elasticity parameter
rho = 1
delta = h/W
gamma = 0.4*delta**2
beta = 1.25             # dimensionless elasticity parameter = lambda/mu
lambda_ = beta          # lame elasticity parameter
g = gamma

N = 64;
# Create mesh and define function space

domain = Rectangle(Point(0.0, 0.0), Point(W, H))
f_domain = Rectangle(Point(0.0, h), Point(W, H))
s_domain = Rectangle(Point(0.0, 0.0), Point(W, h))

domain.set_subdomain(1,s_domain)

mesh = generate_mesh(domain, N)

tol = 1E-14

class Fluid(SubDomain):
    ### Fluid domain is 0 < x < 2.0 and h < y < 2
    def inside(self, x, on_boundary):
        ###return True if 0.0 <= x[0] <= 2.0 and  h <= x[1] <=  2.0 else False
        ##return (between(x[0], (0.0, W)) and between(x[1], (h , H)))
        return True if x[1] >= h - tol else False

class Structure(SubDomain):
    ### Structure domain is 0 < x < 2.0 and 0 < y < h
    def inside(self, x, on_boundary):
        #return True if 0.0 <= x[0] <= 2.0 and 0.0 <= x[1] <=  h else False
        ##return (between(x[0], (0.0, W)) and between(x[1], (0.0, h)))
        return True if x[1] <= h + tol else False

fluid = Fluid()
structure = Structure()

subdomains = CellFunction('size_t', mesh)
fluid.mark(subdomains, 0)
structure.mark(subdomains, 1)

mesh_f = SubMesh(mesh, subdomains, 0)
mesh_s = SubMesh(mesh, subdomains, 1)


#mesh = RectangleMesh(Point(0, 0), Point(W, h), 10, 3)
V = VectorFunctionSpace(mesh_s, 'P', 1)

# Define boundary condition
tol = 1E-14

def clamped_boundary_left(x, on_boundary):
    return on_boundary and x[0] < tol

def clamped_boundary_right(x, on_boundary):
    return on_boundary and x[0] > W - tol

#def clamped_boundary_bottom(x, on_boundary):
    #return on_boundary and x[1] < tol


bc_left = DirichletBC(V, Constant((0, 0)), clamped_boundary_left)

bc_right = DirichletBC(V, Constant((0, 0)), clamped_boundary_right)

bc_bottom = DirichletBC(V, Constant((0, 0)), clamped_boundary_bottom)

bcs = [bc_left, bc_right, bc_bottom]
#bcs = [bc_left, bc_right]


# Mark boundaries

facets = FacetFunction("size_t", mesh_s)
facets.set_all(0)

fsi = CompiledSubDomain('near(x[1], h1)  && on_boundary', h1 = h)
fsi.mark(facets, 3)

dA = Measure('ds', domain = mesh_s, subdomain_data = facets)

# Define strain and stress

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    #return sym(nabla_grad(u))

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
#f = Constant((0, -rho*g))
f = Constant((0, 0))

nodal_sigma = np.loadtxt('nodal_sigma', dtype = float)

T_space = TensorFunctionSpace(mesh_s, 'P', 1)
dofs_s_T = T_space.tabulate_dof_coordinates().reshape((T_space.dim(),-1))
i_s_T = np.where((dofs_s_T[:,1] == h))[0]

sigma_FSI = Function(T_space)
sigma_FSI.vector()[:] = nodal_sigma


dofs_s_V = V.tabulate_dof_coordinates().reshape((V.dim(),-1))
i_s_V_L = np.where((dofs_s_V[:,0] == 0))[0]
i_s_V_R = np.where((dofs_s_V[:,0] == W))[0]
i_s_V_B = np.where((dofs_s_V[:,1] == 0))[0]


# appears as though sigma has values for all of T_space. Isn't really what we want.

# Calculate T...
#n = FacetNormal(mesh_s)
#T = dot(sigma_FSI, n)
#F = I + grad(u)
#T = J*inv(F)*sigma_FSI*N

#T = Constant((-1, 0))      # Traction sigma.n
#T = Expression(('0','x[0]*(x[0]-W)'), W = W, degree = 1)

# create a curve that is 0 at both boundaries and max -ve in middle.

a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*dA(3)

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Plot solution
plot(u, title='Displacement', mode='displacement')

# Plot stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)  # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
V1 = FunctionSpace(mesh_s, 'P', 1)
von_Mises = project(von_Mises, V1)
plot(von_Mises, title='Stress intensity')

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V1)
plot(u_magnitude, 'Displacement magnitude')
print('min/max u:',
      u_magnitude.vector().array().min(),
      u_magnitude.vector().array().max())

# Save solution to file in VTK format
File('elasticity_new/displacement.pvd') << u
File('elasticity_new/von_mises.pvd') << von_Mises
File('elasticity_new/magnitude.pvd') << u_magnitude

# Hold plot
interactive()
