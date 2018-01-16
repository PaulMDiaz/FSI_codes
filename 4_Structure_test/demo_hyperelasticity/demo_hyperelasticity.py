""" This demo program solves a hyperelastic problem. It is implemented
in Python by Johan Hake following the C++ demo by Harish Narayanan"""

# Copyright (C) 2008-2010 Johan Hake and Garth N. Wells
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Harish Narayanan 2009
# Modified by Anders Logg 2011
#
# First added:  2009-10-11
# Last changed: 2012-11-12

# Begin demo

from dolfin import *
import numpy as np

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True}#, \
#               "precompute_ip_const": True}
#               "eliminate_zeros": True, \
#               "precompute_basis_const": True, \

# Create mesh and define function space
mesh = UnitCubeMesh(8, 8, 8)
#mesh = UnitCubeMesh(24, 16, 16)

V_space = VectorFunctionSpace(mesh, "Lagrange", 1)

# Define functions
#du = TrialFunction(V)            # Incremental displacement
#v  = TestFunction(V)             # Test function
#u  = Function(V)                 # Displacement from previous iteration

# Test whether solution can be loaded into mixed function space.
V_element = VectorElement("CG", mesh.ufl_cell() , 1)
mixed_element = FunctionSpace(mesh, MixedElement([V_element, V_element]))

V = TestFunction(mixed_element)
dU = TrialFunction(mixed_element)
U = Function(mixed_element)
U0 = Function(mixed_element)

xi, eta = split(V)
u, v = split(U)

#u0, v0 = U0.split(deepcopy = True)
u0 = Function(V_space)
_u0 = np.loadtxt('twisty_downloaded.txt', dtype = float)
_u0 = np.loadtxt('twisty.txt', dtype = float)

u0.vector()[:] = _u0[:]
#_u02 = np.concatenate((_u0 ,0.0*_u0),axis = 0)

#U.vector()[:] = _u02
#u, v1 = U.split(deepcopy = True)

# Load initial conditions to u0 and v0. Otherwise set to 0.
#u0 = Constant((0,)*V_space.mesh().geometry().dim())
v0 = Constant((0,)*V_space.mesh().geometry().dim())

# Functions for solver
xi, eta = split(V) 	# Test functions
u, v = split(U)		# Functions

cells = CellFunction("size_t", mesh)
dx = Measure('dx', domain = mesh, subdomain_data = cells)

# Project u0 and v0 into U0
a_proj = inner(dU, V)*dx
L_proj = inner(u0, xi)*dx + inner(v0, eta)*dx
solve(a_proj == L_proj, U0)

u0, v0 = U0.split(deepcopy = True)

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# Define Dirichlet boundary (x = 0 or x = 1)
c = Expression(("0.0", "0.0", "0.0"), degree=2)
r = Expression(("scale*0.0",
                "scale*(y0 + (x[1] - y0)*cos(theta) - (x[2] - z0)*sin(theta) - x[1])",
                "scale*(z0 + (x[1] - y0)*sin(theta) + (x[2] - z0)*cos(theta) - x[2])"),
                scale = 0.5, y0 = 0.5, z0 = 0.5, theta = pi/3, degree=2)
#r = Expression(("0.0", "0.0", "0.0"), degree=2)

#clamp = Constant((0.0, 0.0, 0.0))
#bcl = DirichletBC(V, clamp, left)

#bcl = DirichletBC(V, c, left)
#bcr = DirichletBC(V, r, right)

#bcs = [bcl, bcr]
#bcs = [bcl]




B  = Constant((0.0, -0.5, 0.0))  # Body force per unit volume
#B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
#T  = Constant((0.1,  0.1, 0.1))  # Traction force on the boundary
T  = Constant((0.1,  0.0, 0.0))  # Traction force on the boundary

# Kinematics
d = u.geometric_dimension()
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

E = (C - I)/2
E = variable(E)

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Elasticity parameters
#E, nu = 10.0, 0.3
#mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))

mu    = Constant(3.85)
lmbda = Constant(5.77)

# Stored strain energy density (compressible neo-Hookean model)
#psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# (st venant model)
psi = lmbda/2*(tr(E)**2)+mu*tr(E*E)

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Compute first variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Solve variational problem
#solve(F == 0, u, bcs, J=J,
#      form_compiler_parameters=ffc_options)

# Save solution in VTK format
file = File("results/displacement.pvd");

U.assign(U0)
u, v = U.split()
file << u;

# save displacement field
#u_txt = u.vector().array()
#u_txt = u.vector().get_local()

#np.savetxt('twisty.txt', u_txt)

# Plot and hold solution
#plot(u, mode = "displacement", interactive = True)
