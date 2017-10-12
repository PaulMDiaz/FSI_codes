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
nu_s = 0.3              # Poisson's ration
E = 6E6                 # Elastic modulus
rho_s = 1


# first lame constant
mu_s = Constant( E/(2.0*(1.0+nu_s)))
# second lame constant
lambda_s = E*nu_s/((1.0+nu_s)*(1.0-2.0*nu_s))

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
V = VectorFunctionSpace(mesh_s, 'Lagrange', 1)

# Define boundary condition
tol = 1E-14

def clamped_boundary_left(x, on_boundary):
    return on_boundary and x[0] < tol

def clamped_boundary_right(x, on_boundary):
    return on_boundary and x[0] > W - tol

def clamped_boundary_bottom(x, on_boundary):
    return on_boundary and x[1] < tol

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

cells_s = CellFunction('size_t', mesh_s)


dA = Measure('ds', domain = mesh_s, subdomain_data = facets)
dV = Measure('dx', domain = mesh_s, subdomain_data = cells_s)





# Define functions
u = Function(V)             # Displacement from previous iteration
v = TestFunction(V)         # Test function
du = TrialFunction(V)       # Incremental displacement
B = Constant((0.0, 0.0))    # Structure body force.
T = Constant((0.0, 0.0))    # Traction force on the boundary


u0_s = Function(V)
u00_s = Function(V)

Dim = mesh.topology().dim()
I = Identity(Dim)
d_geo = u.geometric_dimension()  # space dimension

# load sigma from fsi problem, project onto structure.
# loaded sigma that has already been projected...
nodal_sigma = np.loadtxt('nodal_sigma', dtype = float)

T_space = TensorFunctionSpace(mesh_s, 'P', 1)
dofs_s_T = T_space.tabulate_dof_coordinates().reshape((T_space.dim(),-1))
i_s_T = np.where((dofs_s_T[:,1] == h))[0]

sigma_FSI = Function(T_space)
sigma_FSI.vector()[:] = nodal_sigma

# should check this...
#sigma_FSI_1 = project(sigma_FSI, T_space, solver_type = "mumps",\
#    form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )

dofs_s_V = V.tabulate_dof_coordinates().reshape((V.dim(),-1))
i_s_V_L = np.where((dofs_s_V[:,0] == 0))[0]
i_s_V_R = np.where((dofs_s_V[:,0] == W))[0]
i_s_V_B = np.where((dofs_s_V[:,1] == 0))[0]

# NeoHookean solver:

# Kinematics
F = I + grad(u) # deformation gradient
C = F.T*F       # Right Cauch-Green tensor

# Invariants of deformation tensor
Ic = tr(C)
J = det(F)

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu_s/2)*(Ic - 3) - mu_s*ln(J) + (lambda_s/2)*(ln(J))**2

# Calculate T...
n = FacetNormal(mesh_s)
#T = dot(sigma_FSI, n)

#F = I + grad(u)
T = J*inv(F)*sigma_FSI*n
#T = dot(T_1, n)

#T = J*sigma_FSI*inv(F).T

#T = Constant((-1, 0))      # Traction sigma.n
#T = Expression(('0','x[0]*(x[0]-W)'), W = W, degree = 1)

# Total potential energy
Pi = psi*dx-dot(T, u)*dA(3) - dot(B, u)*dx
# First directional derivative of Pi about d in the direction of v

Form_s = derivative(Pi, u, v)

# using indices like Abali's code:
i,j, k, l, m = indices(5)
delta = Identity(Dim)
f = Constant((0.0, 0.0))
dt = 0.01

F_s = as_tensor( u[k].dx(i) + delta[k, i], (k, i) )
J_s = det(F_s)

C_s = as_tensor( F_s[k, i]*F_s[k, j], (i, j) )
E_s = as_tensor(1./2.*(C_s[i, j] - delta[i, j]), (i, j) )
S_s = as_tensor( lambda_s*E_s[k, k]*delta[i, j] + 2.*mu_s*E_s[i, j], (i, j))
#S_s = as_tensor( lambda_s*E[k, k]*delta[i, j] + 2.*mu_s*E[i, j], (i, j))
P_s = as_tensor( F_s[i, j]*S_s[k, j], (k, i) )

t_hat = as_tensor( J_s*inv(F_s)[k, j]*sigma_FSI[j, i]*n[k] , (i, ) )


Form_s = ( rho_s*(u-2.*u0_s+u00_s)[i]/(dt*dt)*v[i] + P_s[k, i]*v[i].dx(k) - rho_s*f[i]*v[i] )*dV - \
         t_hat[i]*v[i]*dA(3)

#Form_s = derivative(P_s, u, v)
# Jacobian of the directional derivative Fd
#Gain_s = derivative(Form_s, u, du)

begin("Computing structure displacement")

parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True}

solve(Form_s == 0, u, bcs, J = Gain_s, \
    solver_parameters  = {"newton_solver":{"linear_solver" : "mumps", "relative_tolerance" : 1e-3} }, \
    form_compiler_parameters = {"cpp_optimize" : True, "representation" : "quadrature", "quadrature_degree" : 2} )
#end()


# Plot solution
plot(u, title='Displacement', mode='displacement')

# Plot stress
#s = sigma(d) - (1./3)*tr(sigma(u))*Identity(d_geo)  # deviatoric stress
#von_Mises = sqrt(3./2*inner(s, s))
#V1 = FunctionSpace(mesh_s, 'P', 1)#
#von_Mises = project(von_Mises, V1)
#plot(von_Mises, title='Stress intensity')

V1 = FunctionSpace(mesh_s, 'P', 1)
# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V1)
plot(u_magnitude, 'Displacement magnitude')
print('min/max u:',
      u_magnitude.vector().array().min(),#
      u_magnitude.vector().array().max())

# Save solution to file in VTK format
File('elasticity_neo/displacement.pvd') << u
#File('elasticity_new/von_mises.pvd') << von_Mises
#File('elasticity_new/magnitude.pvd') << u_magnitude

# Hold plot
interactive()
