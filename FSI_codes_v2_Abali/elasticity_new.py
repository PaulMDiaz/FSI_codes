"""
FEniCS tutorial demo program: Linear elastic problem.
  -div(sigma(u)) = f
The model is used to simulate an elastic beam clamped at
its left end and deformed under its own weight.
"""

from __future__ import print_function
from fenics import *

# Scaled variables
W = 2; H = 0.5
mu = 1                  # lame elasticity parameter
rho = 1
delta = H/W
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta          # lame elasticity parameter
g = gamma

# Create mesh and define function space
mesh = RectangleMesh(Point(0, 0), Point(W, H), 10, 3)
V = VectorFunctionSpace(mesh, 'P', 1)

# Define boundary condition
tol = 1E-14

def clamped_boundary_left(x, on_boundary):
    return on_boundary and x[0] < tol

def clamped_boundary_right(x, on_boundary):
    return on_boundary and x[0] > W - tol

bc_left = DirichletBC(V, Constant((0, 0)), clamped_boundary_left)

bc_right = DirichletBC(V, Constant((0, 0)), clamped_boundary_right)

bcs = [bc_left, bc_right]

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
f = Constant((0, -rho*g))
T = Constant((0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Plot solution
plot(u, title='Displacement', mode='displacement')

# Plot stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)  # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)
plot(von_Mises, title='Stress intensity')

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
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
